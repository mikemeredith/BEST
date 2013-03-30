# This file has the S3 generic 'hdi' function and a series of methods.

# hdi.default replaces HDIofMCMC
# hdi.function takes the place of HDIofICDF

hdi <- function(object, credMass=0.95, ...) UseMethod("hdi")

hdi.default <- function(object, credMass=0.95, ...) {
  if(!is.numeric(object))
    stop(paste("No applicable method for class", class(object)))
  if(is.na(credMass) || length(credMass) != 1 || credMass <= 0 || credMass >= 1)
    stop("credMass must be in 0 < credMass < 1")
  if(all(is.na(object)))
    return(c(lower = NA_real_, upper = NA_real_))
  # This is Mike's code from way back:
  x <- sort(object)  # also removes NAs
  n <- length(x)
  # exclude <- ceiling(n * (1 - credMass)) # Not always the same as...
  exclude <- n - floor(n * credMass)       # Number of values to exclude
  low.poss <- x[1:exclude]             # Possible lower limits...
  upp.poss <- x[(n - exclude + 1):n]   # ... and corresponding upper limits
  best <- which.min(upp.poss - low.poss)      # Combination giving the narrowest interval
  result <- c(lower = low.poss[best], upper = upp.poss[best])

  attr(result, "credMass") <- credMass
  return(result)
}

hdi.matrix <- function(object, credMass=0.95, ...) {
  result <- apply(object, 2, hdi.default, credMass=credMass, ...)
  attr(result, "credMass") <- credMass
  return(result)
}

hdi.data.frame <- function(object, credMass=0.95, ...) 
  hdi.matrix(as.matrix(object), credMass=credMass, ...)


hdi.mcmc.list <- function(object, credMass=0.95, ...) 
  hdi.matrix(as.matrix(object), credMass=credMass, ...)

# hdi.bugs <- function(object, credMass=0.95, ...) 
  # hdi.matrix(object$sims.matrix, credMass=credMass, ...)

# hdi.rjags <- function(object, credMass=0.95, ...) 
  # hdi.matrix(object$BUGSoutput$sims.matrix, credMass=credMass, ...)

hdi.function <- function(object, credMass=0.95, tol, ...)  {
  if(is.na(credMass) || length(credMass) != 1 || credMass <= 0 || credMass >= 1)
    stop("credMass must be in 0 < credMass < 1")
  if(missing(tol))
    tol <- 1e-8
  if(class(try(object(0.5, ...), TRUE)) == "try-error")
    stop(paste("Incorrect arguments for the inverse cumulative density function",
        substitute(object)))
  # cf. code in Kruschke 2011 p630
   intervalWidth <- function( lowTailPr , ICDF , credMass , ... ) {
      ICDF( credMass + lowTailPr , ... ) - ICDF( lowTailPr , ... )
   }
   optInfo <- optimize( intervalWidth , c( 0 , 1.0 - credMass) , ICDF=object ,
                        credMass=credMass , tol=tol , ... )
   HDIlowTailPr <- optInfo$minimum
   result <- c(lower = object( HDIlowTailPr , ... ) ,
	            upper = object( credMass + HDIlowTailPr , ... ) )
  attr(result, "credMass") <- credMass
  return(result)
}




