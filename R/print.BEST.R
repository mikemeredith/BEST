print.BEST <- function(x, digits=4, ...) {
  # Somewhat less quick and dirty print method for BEST objects.
  mat <- as.matrix(x)

  # If we decide that BESTmcmc should return a matrix, this could have attributes
  #   Rhat and n.eff.

  Rhat <- attr(x, "Rhat")
  if(is.null(Rhat) && inherits(x, "mcmc.list")) {
    Rhat <- gelman.diag(x)$psrf[, 1]
    attr(x, "Rhat") <- Rhat
  }
  n.eff <- attr(x, "n.eff")
  if(is.null(n.eff) && inherits(x, "mcmc.list")){
    n.eff <- effectiveSize(x)
    attr(x, "n.eff") <- n.eff
  }

  toPrint <- cbind(
    mean = colMeans(mat),
    sd = apply(mat, 2, sd),
    median = apply(mat, 2, median), 
    t(hdi(mat)))
  colnames(toPrint)[4:5] <- c("HDIlo", "HDIup")
  if(!is.null(Rhat))
    toPrint <- cbind(toPrint, Rhat = Rhat)
  if(!is.null(n.eff))
    toPrint <- cbind(toPrint, n.eff = round(n.eff))

  cat("MCMC fit results for BEST analysis:\n")
  cat(nrow(mat), "simulations saved.\n")
  print(toPrint, digits = digits)
  cat("\n'HDIlo' and 'HDIup' are the limits of a 95% HDI credible interval.\n")
  if(!is.null(Rhat))
    cat("'Rhat' is the potential scale reduction factor (at convergence, Rhat=1).\n")
  if(!is.null(n.eff))
    cat("'n.eff' is a crude measure of effective sample size.\n")

}
