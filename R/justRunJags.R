
# Functions to run JAGS via 'rjags' with no extra features/annoyances

# data, params, modelFile have the usual meanings
# initList must be a list with one component per chain (NOT a function)
# sample = the number of values required
# adapt, burnin - note that adaptation contines during the burn-in phase.

# Run JAGS in serial mode.

# This function is also called (with chains=1) to run JAGS in each worker.
# Note that initList MUST be the first argument to work with parLapply.
justRunJagsSerial <- function(initList, data, params, modelFile,
    chains, sample, burnin, adapt=1000, thin=1) {
  jm <- rjags::jags.model(modelFile, data, initList, n.chains=chains, n.adapt=0)
  update(jm, adapt + burnin)
  if(!rjags::adapt(jm, 0, end.adaptation=TRUE))
    warning("Adaptation was not adequate.")
  cat("\nSampling from the posterior distributions:\n")
  rjags::coda.samples(jm, params, n.iter=ceiling(sample / chains) * thin, thin=thin)
}
# ---------------------------------------------------------------

# The main function to run JAGS

justRunJags <- function(data, initList, params, modelFile,
        chains, sample, burnin, thin=1, adapt = 1000,
        parallel = NULL, seed=NULL, verbose=verbose)  {

  if(parallel) {   ##### Do the parallel stuff #####
    if(verbose) {
      message("Waiting for parallel processing to complete...", appendLF=FALSE)
      flush.console()
    }
    cl <- makeCluster(3) ; on.exit(stopCluster(cl))
    clusterEvalQ(cl, library(rjags))
    chainList <- parLapply(cl, initList, justRunJagsSerial, data=data, params=params,
      modelFile=modelFile, chains=1, sample=ceiling(sample / chains), burnin=burnin, adapt=adapt, thin=thin)
    mcmcList <- coda::mcmc.list(lapply(chainList, function(x) x[[1]]))
    if(verbose)
      message("done.")
  } else {     ##### Do the serial stuff #####
    if(verbose) {
      mcmcList <- justRunJagsSerial(initList, data=data, params=params, modelFile=modelFile,
                  chains=chains, sample=sample, burnin=burnin, adapt=adapt, thin=thin)
    } else {
      null <- capture.output(
        mcmcList <- justRunJagsSerial(initList, data=data, params=params, modelFile=modelFile,
                  chains=chains, sample=sample, burnin=burnin, adapt=adapt, thin=thin) )
    }
  }

  invisible(mcmcList)
}
