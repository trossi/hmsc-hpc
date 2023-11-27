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

#### Step 2. Export initial model ####

set.seed(RS+42)
init_obj = sampleMcmc(m, samples=nSamples, thin=thin,
                      transient=transient,
                      nChains=nChains, verbose=verbose, engine="pass")

init_file_name = sprintf("TF-init-obj-%s.rds", selected_experiment$name)
init_file_path = file.path(path, "examples/data", init_file_name)
python_file_name = "run_gibbs_sampler.py"
python_file_path = file.path(path, "examples", python_file_name)
postList_file_name = sprintf("TF-postList-obj-%s.rds", selected_experiment$name)
postList_file_path = file.path(path, "examples/data", postList_file_name)

nr = init_obj[["hM"]][["nr"]]
rLNames = init_obj[["hM"]][["ranLevelsUsed"]]
for (r in seq_len(nr)) {
  rLName = rLNames[[r]]
  init_obj[["hM"]][["rL"]][[rLName]][["s"]] = NULL
  init_obj[["hM"]][["ranLevels"]][[rLName]][["s"]] = NULL
  spatialMethod = init_obj[["hM"]][["rL"]][[r]][["spatialMethod"]]
  if (!is.null(spatialMethod)) {
    if (spatialMethod == "NNGP") {
      gN = length(init_obj[["dataParList"]][["rLPar"]][[r]][["iWg"]])
  
      for (i in seq_len(gN)) {
        iWg = as(init_obj[["dataParList"]][["rLPar"]][[r]][["iWg"]][[i]], "dgTMatrix")
        RiWg = as(init_obj[["dataParList"]][["rLPar"]][[r]][["RiWg"]][[i]], "dgTMatrix")
        init_obj[["dataParList"]][["rLPar"]][[r]][["iWgi"]][[i]] = iWg@i
        init_obj[["dataParList"]][["rLPar"]][[r]][["iWgj"]][[i]] = iWg@j
        init_obj[["dataParList"]][["rLPar"]][[r]][["iWgx"]][[i]] = iWg@x
        init_obj[["dataParList"]][["rLPar"]][[r]][["RiWgi"]][[i]] = RiWg@i
        init_obj[["dataParList"]][["rLPar"]][[r]][["RiWgj"]][[i]] = RiWg@j
        init_obj[["dataParList"]][["rLPar"]][[r]][["RiWgx"]][[i]] = RiWg@x
      }
      init_obj[["dataParList"]][["rLPar"]][[r]][["iWg"]] = NULL
      init_obj[["dataParList"]][["rLPar"]][[r]][["RiWg"]] = NULL
    }
    else if (spatialMethod == "GPP") {
      init_obj[["dataParList"]][["rLPar"]][[r]][["nK"]] = nrow(init_obj[["dataParList"]][["rLPar"]][[1]][["Fg"]])
    }
  }
}

# Replace hM object with json representation for export
hM = init_obj$hM
init_obj$hM <- to_json(hM)

saveRDS(init_obj, file = init_file_path, compress=TRUE)

# Restore hM object
init_obj$hM <- hM


#### Step 4. Run TF code ####

python_cmd = paste("python", sprintf("'%s'",python_file_path),
                   "--samples", nSamples,
                   "--transient", transient,
                   "--thin", thin,
                   "--verbose", verbose,
                   "--input", init_file_path,
                   "--output", postList_file_path)
print(python_cmd)
