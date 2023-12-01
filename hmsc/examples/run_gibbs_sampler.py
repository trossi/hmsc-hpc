import os
from contextlib import nullcontext
import time
import argparse
import numpy as np
import tensorflow as tf
from hmsc.gibbs_sampler import GibbsSampler
from hmsc.utils.export_rds_utils import (
    load_model_from_rds,
    save_chains_postList_to_rds,
)
from hmsc.utils.import_utils import (
    load_model_dims,
    load_model_data,
    load_prior_hyperparams,
    load_random_level_hyperparams,
    load_model_hyperparams,
    init_params,
)
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"


def load_params(file_path, dtype=np.float64):
    hmscImport, hmscModel = load_model_from_rds(file_path)
    modelDims = load_model_dims(hmscModel)
    modelData = load_model_data(hmscModel, hmscImport.get("initParList"))
    priorHyperparams = load_prior_hyperparams(hmscModel)
    # currently not used at all
    # modelHyperparams = load_model_hyperparams(hmscModel, hmscImport.get("dataParList"))
    modelHyperparams = None
    rLHyperparams = load_random_level_hyperparams(hmscModel, hmscImport.get("dataParList"))
    initParList = init_params(hmscImport.get("initParList"), modelData, modelDims, rLHyperparams)
    nChains = int(hmscImport.get("nChains")[0])
  
    return modelDims, modelData, priorHyperparams, modelHyperparams, rLHyperparams, initParList, nChains


def run_gibbs_sampler(
    num_samples,
    sample_thining,
    sample_burnin,
    verbose,
    init_obj_file_path,
    postList_file_path,
    chainIndList=None,
    truncated_normal_library="tf",
    flag_save_eta=True,
    flag_save_postList_to_rds=True,
    flag_profile=False,
):
    (
        modelDims,
        modelData,
        priorHyperparams,
        modelHyperparams, #this precomputed one (e.g. Qg) is currently not used and is computed at runtime
        rLHyperparams,
        initParList,
        nChainsTotal,
    ) = load_params(init_obj_file_path)
    gibbs = GibbsSampler(modelDims, modelData, priorHyperparams, rLHyperparams)
    
    if chainIndList is None:
      chainIndList = [*range(nChainsTotal)]
    else:
      chainIndListNew = [chainInd for chainInd in chainIndList if chainInd < nChainsTotal]
      if chainIndList != chainIndListNew:
        print("Input chainIndList", chainIndList)
        print("Adjusted chainIndList", chainIndListNew)
        chainIndList = chainIndListNew

    print("Initializing TF graph", flush=True)
    startTime = time.time()
    parSamples = gibbs.sampling_routine(
        initParList[0],
        num_samples=tf.constant(1),
        sample_burnin=tf.constant(1),
        sample_thining=tf.constant(1),
        verbose=verbose,
        truncated_normal_library=truncated_normal_library,
        flag_save_eta=flag_save_eta,
    )
    elapsedTime = time.time() - startTime
    print("TF graph initialized in %.1f sec" % elapsedTime, flush=True)
     
    print("Running TF Gibbs sampler for %d chains with indices" % len(chainIndList), chainIndList, flush=True)
    with tf.profiler.experimental.Profile('logdir') if flag_profile else nullcontext():
        startTime = time.time()
        postList = [None] * len(chainIndList)
        
        for chainInd, chain in enumerate(chainIndList):
            print("\n", "Computing chain %d" % chain, flush=True)
    
            parSamples = gibbs.sampling_routine(
                initParList[chain],
                num_samples=tf.constant(num_samples),
                sample_burnin=tf.constant(sample_burnin),
                sample_thining=tf.constant(sample_thining),
                verbose=verbose,
                truncated_normal_library=truncated_normal_library,
                flag_save_eta=flag_save_eta,
            )
            postList[chainInd] = [None] * num_samples
            for n in range(num_samples):
                parSnapshot = {
                    "Beta": parSamples["Beta"][n],
                    "BetaSel": [samples[n] for samples in parSamples["BetaSel"]],
                    "Gamma": parSamples["Gamma"][n],
                    "iV": parSamples["iV"][n],
                    "rhoInd": parSamples["rhoInd"][n],
                    "sigma": parSamples["sigma"][n],
                    "Lambda": [samples[n] for samples in parSamples["Lambda"]],
                    "Psi": [samples[n] for samples in parSamples["Psi"]],
                    "Delta": [samples[n] for samples in parSamples["Delta"]],
                    "Eta": [samples[n] for samples in parSamples["Eta"]] if flag_save_eta else None,
                    "AlphaInd": [samples[n] for samples in parSamples["AlphaInd"]],
                    "wRRR": parSamples["wRRR"][n] if "wRRR" in parSamples else None,
                    "PsiRRR": parSamples["PsiRRR"][n] if "PsiRRR" in parSamples else None,
                    "DeltaRRR": parSamples["DeltaRRR"][n] if "DeltaRRR" in parSamples else None,
                }
                postList[chainInd][n] = parSnapshot
            
            elapsedTime = time.time() - startTime
            print("\n", "%d chains completed in %.1f sec" % (chainInd+1, elapsedTime), flush=True)
    
        elapsedTime = time.time() - startTime
        print("\n", "Whole Gibbs sampler elapsed %.1f" % elapsedTime, flush=True)
    
    if flag_save_postList_to_rds:
        save_chains_postList_to_rds(postList, postList_file_path, len(chainIndList), elapsedTime, flag_save_eta)


if __name__ == "__main__":

    argParser = argparse.ArgumentParser()
    argParser.add_argument(
        "-s",
        "--samples",
        type=int,
        default=2,
        help="number of samples obtained per chain",
    )
    argParser.add_argument(
        "-b",
        "--transient",
        type=int,
        default=0,
        help="number of samples discarded before recording posterior samples",
    )
    argParser.add_argument(
        "-t",
        "--thin",
        type=int,
        default=1,
        help="number of samples between each recording of posterior samples",
    )
    argParser.add_argument(
        "-c",
        "--chains",
        type=int,
        nargs='+',
        default=None,
        help="indices of chains to fit",
    )
    argParser.add_argument(
        "-i",
        "--input",
        type=str,
        default="TF-init-obj.rds",
        help="input RDS file with parameters for model initialization",
    )
    argParser.add_argument(
        "-o",
        "--output",
        type=str,
        default="TF-postList-obj.rds",
        help="output RDS file with recorded posterier samples",
    )
    argParser.add_argument(
        "-v",
        "--verbose",
        type=int,
        default=1,
        help="print out information meassages and progress status",
    )
    argParser.add_argument(
        "--tnlib",
        type=str,
        default="tf",
        choices=["scipy", "tf", "tfd"],
        help="which library is used for sampling trunacted normal: scipy, tf or tfd",
    )
    argParser.add_argument(
        "--fse",
        type=int,
        default=1,
        help="whether to save Eta posterior",
    )
    argParser.add_argument(
        "--profile",
        type=int,
        default=0,
        help="whether to run profiler alongside sampling",
    )
    args = argParser.parse_args()
    print("args=%s" % args)
    print("working directory", os.getcwd(), flush=True)
    init_obj_file_path = args.input
    postList_file_path = args.output

    run_gibbs_sampler(
        num_samples=args.samples,
        sample_thining=args.thin,
        sample_burnin=args.transient,
        verbose=args.verbose,
        init_obj_file_path=init_obj_file_path,
        postList_file_path=postList_file_path,
        chainIndList=args.chains,
        truncated_normal_library=args.tnlib,
        flag_save_eta=bool(args.fse),
        flag_save_postList_to_rds=True,
        flag_profile=bool(args.profile),
    )
