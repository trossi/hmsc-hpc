library(Hmsc)
library(jsonify)
library(vioplot)
library(abind)
RS = 1
set.seed(RS)
nParallel = 8
flagFitR = 0

path = getwd()

#### Step 1. Load model ####

load(file = file.path(path, "examples/data", "unfitted_models_2.RData"))
experiments = list(
  M1=list(name="model1",id=1),
  M2=list(name="model2",id=2),
  M3=list(name="model3",id=3),
  M4=list(name="model4",id=4),
  M5=list(name="model5",id=5),
  M6=list(name="model6",id=6),
  M7=list(name="model7",id=7),
  M8=list(name="model8",id=8),
  M9=list(name="model9",id=9),
  M10=list(name="model10",id=10),
  M11=list(name="model11",id=NA)
)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
    stop("Missing args")
}
print(args)
selected_experiment = experiments[[sprintf("M%s", args[1])]]
if (is.null(selected_experiment)) {
    stop("Unknown experiment")
}

if(!is.na(selected_experiment$id)){
  m = models[[selected_experiment$id]]
}
print(m)

# nChains = 8
# nSamples = 250
# thin = 100
nChains = 8
nSamples = 3
thin = 5

transient = nSamples*thin
verbose = thin*1


#### Step 3. Run R code ####

#
# Generate sampled posteriors in R
#

if(flagFitR){
  set.seed(RS+42)
  startTime = proc.time()
  obj.R = sampleMcmc(m, samples = nSamples, thin = thin,
                     transient = transient, 
                     nChains = nChains, nParallel=min(nChains,nParallel),
                     verbose = verbose, updater=list(Gamma2=FALSE, GammaEta=FALSE)) #fitted by R
  elapsedTime = proc.time() - startTime
  print(elapsedTime)
  save(obj.R, elapsedTime, file=fitR_file_path)
} else{
  load(file=fitR_file_path)
}


#### Step 5. Import TF posteriors ####################################################################################

postList.TF <- from_json(readRDS(file = postList_file_path)[[1]])

names(postList.TF) = NULL
for (chain in seq_len(nChains)) {
  names(postList.TF[[chain]]) = NULL
}

obj.TF = init_obj$hM
obj.TF[["postList"]] = postList.TF

# tmp = obj.R$postList[[1]]
# obj.R$postList[[1]] = obj.R$postList[[2]]
# tmp = obj.TF$postList[[1]]
# obj.TF$postList[[1]] = obj.TF$postList[[2]]

obj.TF$samples = nSamples
obj.TF$thin = thin
obj.TF$transient = transient

#
# Rescaling Beta/Gamma; copied from combineParameters.R; need to revisit this section
#
nt = obj.TF[["nt"]]
TrInterceptInd = obj.TF[["TrInterceptInd"]]
TrScalePar = obj.TF[["TrScalePar"]]
ncNRRR = obj.TF[["ncNRRR"]]
ncRRR = obj.TF[["ncRRR"]]
ncsel = obj.TF[["ncsel"]]
XScalePar = obj.TF[["XScalePar"]]
XInterceptInd = obj.TF[["XInterceptInd"]]
XRRRScalePar = obj.TF[["XRRRScalePar"]]

for (chain in seq_len(nChains)) {
  for (sample in seq_len(nSamples)) {
    Beta = obj.TF[["postList"]][[chain]][[sample]][["Beta"]]
    BetaSel = obj.TF[["postList"]][[chain]][[sample]][["BetaSel"]]
    if(is.matrix(BetaSel)){
      BetaSel = split(BetaSel, rep(1:nrow(BetaSel), ncol(BetaSel)))
    }
    Gamma = obj.TF[["postList"]][[chain]][[sample]][["Gamma"]]
    iV = obj.TF[["postList"]][[chain]][[sample]][["iV"]]
    rho = obj.TF$rhopw[obj.TF[["postList"]][[chain]][[sample]][["rhoInd"]], 1]
    sigma = obj.TF[["postList"]][[chain]][[sample]][["sigma"]]

    for(p in 1:nt){
      me = TrScalePar[1,p]
      sc = TrScalePar[2,p]
      if(me!=0 || sc!=1){
        Gamma[,p] = Gamma[,p]/sc
        if(!is.null(TrInterceptInd)){
          Gamma[,TrInterceptInd] = Gamma[,TrInterceptInd] - me*Gamma[,p]
        }
      }
    }

    for(k in 1:ncNRRR){
      me = XScalePar[1,k]
      sc = XScalePar[2,k]
      if(me!=0 || sc!=1){
        Beta[k,] = Beta[k,]/sc
        Gamma[k,] = Gamma[k,]/sc
        if(!is.null(XInterceptInd)){
          Beta[XInterceptInd,] = Beta[XInterceptInd,] - me*Beta[k,]
          Gamma[XInterceptInd,] = Gamma[XInterceptInd,] - me*Gamma[k,]
        }
        iV[k,] = iV[k,]*sc
        iV[,k] = iV[,k]*sc
      }
    }

    for(k in seq_len(ncRRR)){
      me = XRRRScalePar[1,k]
      sc = XRRRScalePar[2,k]
      if(me!=0 || sc!=1){
        Beta[ncNRRR+k,] = Beta[ncNRRR+k,]/sc
        Gamma[ncNRRR+k,] = Gamma[ncNRRR+k,]/sc
        if(!is.null(XInterceptInd)){
          Beta[XInterceptInd,] = Beta[XInterceptInd,] - me*Beta[ncNRRR+k,]
          Gamma[XInterceptInd,] = Gamma[XInterceptInd,] - me*Gamma[ncNRRR+k,]
        }
        iV[ncNRRR+k,] = iV[ncNRRR+k,]*sc
        iV[,ncNRRR+k] = iV[,ncNRRR+k]*sc
      }
    }

    for (i in seq_len(ncsel)){
      XSel = obj.TF$XSelect[[i]]
      for (spg in 1:length(XSel$q)){
        if(!BetaSel[[i]][spg]){
          fsp = which(XSel$spGroup==spg)
          Beta[XSel$covGroup,fsp]=0
        }
      }
    }

    obj.TF[["postList"]][[chain]][[sample]][["Beta"]] = Beta
    obj.TF[["postList"]][[chain]][[sample]][["Gamma"]] = Gamma
    obj.TF[["postList"]][[chain]][[sample]][["V"]] = chol2inv(chol(iV))
    obj.TF[["postList"]][[chain]][[sample]][["iV"]] = NULL
    obj.TF[["postList"]][[chain]][[sample]][["rho"]] = rho
    obj.TF[["postList"]][[chain]][[sample]][["sigma"]] = sigma^2
  }
}
obj.TF = alignPosterior(obj.TF)

#### Step 6. Visualize ####

#
# Plot sampled posteriors summaries from R only and TF
#
# obj.TF = obj.R #CHECK THIS

b.R = getPostEstimate(obj.R, "Beta")
b.TF = getPostEstimate(obj.TF, "Beta")
par(mfrow=c(ceiling(nc/2),min(2,nc-1)))
for(k in 1:obj.R$nc){
  plot(b.R$mean[k,], b.TF$mean[k,], main=paste(k, obj.R$covNames[k]))
  abline(0,1)
  # text(b.R$mean[k,], b.TF$mean[k,])
}
par(mfrow=c(1,1))
if(any(obj.R$distr[,2] != 0)){
  s.R = getPostEstimate(obj.R, "sigma")
  s.TF = getPostEstimate(obj.TF, "sigma")
  plot(s.R$mean, s.TF$mean, main="sigma")
  abline(0,1)
}
if(obj.R$nr > 0){
  omega.R = getPostEstimate(obj.R, "Omega")$mean
  omega.TF = getPostEstimate(obj.TF, "Omega")$mean
  # plot(crossprod(Lambda), omega.R, ylim=range(c(omega.R,omega.TF)), main="Omega")
  # points(crossprod(Lambda), omega.TF, col="blue")
  # plot(crossprod(Lambda)[-19,-19], omega.R[-19,-19], ylim=range(c(omega.R[-19,-19],omega.TF[-19,-19])), main="Omega")
  # points(crossprod(Lambda)[-19,-19], omega.TF[-19,-19], col="blue")
  plot(omega.R, omega.TF, main="Omega")
  abline(0,1)
  
  lambda.R = getPostEstimate(obj.R, "Lambda")$mean[1:2,]
  lambda.TF = getPostEstimate(obj.TF, "Lambda")$mean[1:2,]
  # plot(lambda.R[1,], lambda.TF[2,])
  # plot(lambda.R[2,], lambda.TF[1,])
}


obj.R.TF = obj.R
obj.R.TF$postList = c(obj.R$postList,obj.TF$postList)
obj.list = list(R=obj.R,TF=obj.TF,R.TF=obj.R.TF)

#beta, gamma, sigma
varVec = c("Beta","Gamma")
if(any(obj.R$distr[,2] != 0)) varVec = c(varVec, "Sigma")
for(variable in 1:length(varVec)){
  for(i in 1:3){
    mpost = convertToCodaObject(obj.list[[i]], Lambda=FALSE, Omega=FALSE, Psi=FALSE, Delta=FALSE, Eta=FALSE)
    mpost.var = mpost[[varVec[variable]]]
    psrf = gelman.diag(mpost.var, multivariate=FALSE)$psrf
    if(i == 1) {ma = psrf[,1]} else {ma = cbind(ma,psrf[,1])}
  }
  par(mfrow=c(2,1))
  print(paste(varVec[variable], "NaN", paste(apply(is.nan(ma), 2, sum),collapse=" "), sep=" - "))
  ma[is.nan(ma)] = 0.9
  vioplot(ma,names=names(obj.list),ylim=c(0.9,max(ma, na.rm=TRUE)),main=varVec[variable])
  vioplot(ma,names=names(obj.list),ylim=c(0.9,1.1),main=varVec[variable])
  par(mfrow=c(1,1))
}

#rho
if(!is.null(obj.R$C)){
  mat = rbind(unlist(getPostEstimate(obj.R, "rho")), unlist(getPostEstimate(obj.TF, "rho")))
  rownames(mat) = c("R", "TF")
  cat("rho posterior summary\n")
  print(mat)
}

#omega
maxOmega = 1000 #number of species pairs to be subsampled
nr = obj.list[[1]]$nr
for(k in seq_len(nr)){
  for(i in 1:3){
    mpost = convertToCodaObject(obj.list[[i]], Beta=FALSE, Lambda=FALSE, Psi=FALSE, Delta=FALSE, Eta=FALSE)
    tmp = mpost$Omega[[k]]
    z = dim(tmp[[1]])[2]
    if(z > maxOmega){
      if(i==1) sel = sample(1:z, size = maxOmega)
      for(j in 1:length(tmp)){
        tmp[[j]] = tmp[[j]][,sel]
      }
    }
    psrf = gelman.diag(tmp, multivariate = FALSE)$psrf
    if(i == 1) {ma = psrf[,1]} else {ma = cbind(ma,psrf[,1])}
  }
  par(mfrow=c(2,1))
  print(paste(varVec[variable], sprintf("Omega%d",k), paste(apply(is.nan(ma), 2, sum),collapse=" "), sep=" - "))
  print(paste(varVec[variable], sprintf("Omega%d",k), paste(apply(is.infinite(ma), 2, sum),collapse=" "), sep=" - "))
  ma[is.nan(ma) | is.infinite(ma)] = 0.9
  vioplot(ma,names=names(obj.list),ylim=c(0.9,max(ma)),main=paste("omega",names(obj.list[[1]]$ranLevels)[k]))
  vioplot(ma,names=names(obj.list),ylim=c(0.9,1.1),main=paste("omega",names(obj.list[[1]]$ranLevels)[k]))
  par(mfrow=c(1,1))
}

#reduced rank regression - effective beta
if(obj.R$ncRRR > 0){
  for(i in 1:3){
    BetaRRR_list = vector("list", nChains)
    for (chain in seq_len(nChains)) {
      BetaRRR_samples = vector("list", nSamples)
      for (sample in seq_len(nSamples)) {
        Beta = obj.list[[i]][["postList"]][[chain]][[sample]][["Beta"]]
        w = obj.list[[i]][["postList"]][[chain]][[sample]][["wRRR"]]
        BetaRRR_samples[[sample]] = as.vector(crossprod(w, Beta[-(1:obj.R$ncNRRR),]))
      }
      BetaRRR_list[[chain]] = mcmc(abind(BetaRRR_samples, along=0))
    }
    BetaRRR_mcmc.list = as.mcmc.list(BetaRRR_list)
    meanVal = colMeans(as.matrix(BetaRRR_mcmc.list))
    psrf = gelman.diag(BetaRRR_mcmc.list, multivariate=FALSE)$psrf
    if(i == 1) {ma = psrf[,1]; mv = meanVal} else {ma = cbind(ma,psrf[,1]); mv = cbind(mv,meanVal)}
  }
  plot(mv[,1], mv[,2])
  par(mfrow=c(2,1))
  vioplot(ma,names=names(obj.list),ylim=c(0.9,max(ma)),main="BetaRRR")
  vioplot(ma,names=names(obj.list),ylim=c(0.9,1.1),main="BetaRRR")
  par(mfrow=c(1,1))
}

mpost_R = convertToCodaObject(obj.R, Lambda=FALSE, Omega=FALSE, Psi=FALSE, Delta=FALSE, Eta=FALSE)$V
plot(mpost_R)
mpost_TF = convertToCodaObject(obj.TF, Lambda=FALSE, Omega=FALSE, Psi=FALSE, Delta=FALSE, Eta=FALSE)$V
plot(mpost_TF)
par(mfrow=c(1,1))

# psrf_R = gelman.diag(mpost_R, multivariate=FALSE)$psrf
# psrf_TF = gelman.diag(mpost_TF, multivariate=FALSE)$psrf
# plot(psrf_R[,1], psrf_TF[,1], type="n")
# text(psrf_R[,1], psrf_TF[,1], 1:nrow(psrf_R))



# library(truncnorm)
# for(cInd in 1:chain){
#   # ind = order(rowMeans(abs(obj.R$postList[[cInd]][[nSamples]]$Lambda[[1]])^2), decreasing=TRUE)[1:nc]
#   cInd2 = 1
#   ind = which(rowSums(obj.R$postList[[cInd2]][[nSamples]]$Psi[[1]]) > 0)
#   init_obj$initParList[[cInd]]$Eta[[1]] = obj.R$postList[[cInd2]][[nSamples]]$Eta[[1]][,ind,drop=FALSE]
#   init_obj$initParList[[cInd]]$Lambda[[1]] = obj.R$postList[[cInd2]][[nSamples]]$Lambda[[1]][ind,,drop=FALSE]
#   init_obj$initParList[[cInd]]$Delta[[1]] = 1+0*obj.R$postList[[cInd2]][[nSamples]]$Delta[[1]][ind,,drop=FALSE]
#   init_obj$initParList[[cInd]]$Psi[[1]] = 1+0*obj.R$postList[[cInd2]][[nSamples]]$Psi[[1]][ind,,drop=FALSE]
#   init_obj$initParList[[cInd]]$Alpha[[1]] = obj.R$postList[[cInd2]][[nSamples]]$Alpha[[1]][ind,drop=FALSE]
#   L = init_obj$initParList[[cInd]]$Eta[[1]] %*% init_obj$initParList[[cInd]]$Lambda[[1]]
#   lB = rep(-Inf, length(Y))
#   uB = rep(Inf, length(Y))
#   lB[Y] = 0
#   uB[!Y] = 0
#   z = rtruncnorm(length(Y), a=lB, b=uB, mean=L, sd=1)
#   init_obj$initParList[[cInd]]$Z = matrix(z, nrow(L), ncol(L))
# }
# write(to_json(init_obj), file = init_file_path)


