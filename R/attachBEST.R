# File contains attachBEST and detachBEST.
# (detachBEST moved here 2013-02-15)

attachBEST <-
function(BESTobj, overwrite = NA) 
{
  mcmcChain <- as.matrix(BESTobj)
  # Sanity checks:
  if(!(ncol(mcmcChain) == 5 && all(colnames(mcmcChain) == c("mu[1]", "mu[2]","nu","sigma[1]","sigma[2]"))) &&
     !(ncol(mcmcChain) == 3 && all(colnames(mcmcChain) == c("mu","nu","sigma"))) )
        stop("BESTobj is not a valid BEST object")
  # Turn mcmcChain into a data frame
  BESTsims <- as.data.frame(mcmcChain)
  if(ncol(BESTsims) == 5)
    colnames(BESTsims) <- c("mu1", "mu2", "nu", "sigma1", "sigma2")

  # Remove masking objects in .GlobalEnv
  rem <- colnames(BESTsims)[names(BESTsims) %in% ls(.GlobalEnv)]
  if (length(rem) == 0) {
      overwrite <- FALSE
  } else if (is.na(overwrite)) {
      if (interactive()) {
        question <- paste("The following objects in .GlobalEnv will mask\nobjects in the attached database:\n", 
          paste(rem, collapse = ", "), "\nRemove these objects from .GlobalEnv?", 
          sep = "")
        if (.Platform$OS.type == "windows") {
          overwrite <- "YES" == winDialog(type = "yesno", 
            question)
        } else {
          overwrite <- 1 == menu(c("YES", "NO"), graphics = FALSE, 
            title = question)
        }
      } else {
        overwrite <- FALSE
      }
  }
  if (overwrite) 
    remove(list = rem, envir = .GlobalEnv)
  # Detach any existing BESTsims object:
  while("BESTsims" %in% search())
    detach("BESTsims")

  # Now attach BESTsims
  cat("You can now access", paste(colnames(BESTsims), collapse=", "), "by name.\n")
  attach(BESTsims)
}

detachBEST <-
function() 
{
  while("BESTsims" %in% search()) {
    detach("BESTsims")
  }
}

