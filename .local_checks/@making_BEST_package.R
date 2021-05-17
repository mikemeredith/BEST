
setwd("D:/Github/BEST_package") # my desktop
setwd("C:/Github/BEST_package") # my laptop
dir()

library(spelling)
update_wordlist(pkg = "BEST", confirm = TRUE)
out <- spell_check_package(pkg = "BEST")

devtools::load_all("C:/GitHub/BEST_package/BEST")
system("R CMD INSTALL BEST") # Use this for a "dev" install.

# To check the current CRAN version
# ---------------------------------
# download.packages("BEST", destdir=".", type="source")
# pkg <- "BEST_0.5.2.tar.gz"

# Create the BEST package
# ==========================
unlink(list.files(pattern="Rplots.pdf", recursive=TRUE))
system("R CMD build BEST")  # Produces the .tar.gz
pkg <- "BEST_0.5.3.9000.tar.gz"  # <-- fix version number here ################

# Pick one to check:
## on desktop
system(paste("R CMD check ", pkg))
system(paste("R CMD check ", pkg, "--as-cran"))  # as-cran now runs donttest
## on laptop
system(paste("R CMD check ", pkg, "--no-manual"))
system(paste("R CMD check ", pkg, "--as-cran --no-manual"))

# Pick one to install
system(paste("R CMD INSTALL ", pkg))            # install only
system(paste("R CMD INSTALL ", pkg, "--build")) # install and produce the .zip binary



library(testthat)
test_package("BEST", reporter=ProgressReporter)

# Try it out:
library(BEST)
?BEST

# Run these examples, we need the output:
example("BEST-package")
example(BESTpower)

# Check that the power plots come up ok and results are saved:
unlink("testSave.Rda")  # delete any old testSave.Rda files
system.time(
BESTpower(BESTout, N1=length(y1), N2=length(y2),
            ROPEm=c(-0.1,0.1), maxHDIWm=2.0, nRep=2,
            saveName = "testSave.Rda", showFirst=2, verbose=1)
)
load("testSave.Rda")
power
graphics.off() # Clean up
unlink("testSave.Rda")  # Clean up

# Check that the plots look ok (auto-checks can't do that):
example(plotPost)
example(plotAreaInROPE)

# Check parallel and verbose options
y1 <- c(5.77, 5.33, 4.59, 4.33, 3.66, 4.48)
y2 <- c(3.88, 3.55, 3.29, 2.59, 2.33, 3.59)

system.time(
  BESToutP0 <- BESTmcmc(y1, y2, rnd.seed=123) )  # default (depends on machine)
system.time(
  BESToutPT <- BESTmcmc(y1, y2, rnd.seed=123, parallel=TRUE) )  # 5 secs
system.time(
  BESToutPF <- BESTmcmc(y1, y2, rnd.seed=123, parallel=FALSE) )  # 10 secs
BESToutP0
BESToutPT
BESToutPF

( BESToutP0Q <- BESTmcmc(y1, y2, rnd.seed=123, verbose=FALSE) )
( BESToutPQ <- BESTmcmc(y1, y2, rnd.seed=123, verbose=FALSE, parallel=TRUE) )
( BESToutQ <- BESTmcmc(y1, y2, rnd.seed=123, verbose=FALSE, parallel=FALSE) )

# Check priors only
( BESToutPrior <- BESTmcmc(y1, y2, priors=list(), doPriorsOnly=TRUE, rnd.seed=123) )
plotAll(BESToutPrior)
priors <- list(muSD = 10, sigmaMode = 10, sigmaSD = 100)
( BESToutPriorInf <- BESTmcmc(y1, y2, priors=priors, doPriorsOnly=TRUE, rnd.seed=123) )
plotAll(BESToutPriorInf)


# Check small samples
y1 <- 4
y2 <- 5
( BESToutPs <- BESTmcmc(y1, y2, rnd.seed=123) )
( BESTouts <- BESTmcmc(y1, y2, rnd.seed=123, parallel=FALSE) )
plot(BESToutPs)
plotPostPred(BESToutPs)
plotAll(BESToutPs)
priors <- list(muSD = 100, sigmaMode = 10, sigmaSD = 100)
( BESTout2pi <- BESTmcmc(y1, y2, priors=priors, rnd.seed=123) )
plotAll(BESTout2pi)

y <- 4:5
( BESTout1 <- BESTmcmc(y, rnd.seed=123) )
( BESTout1p <- BESTmcmc(y, priors=list(), rnd.seed=123) )
plot(BESTout1p)
plotPostPred(BESTout1p)
priors <- list(muSD = 10, sigmaMode = 10, sigmaSD = 100)
( BESTout1pi <- BESTmcmc(y, priors=priors, rnd.seed=123) )
plotAll(BESTout1pi)
y <- 4
( BESTout1pi <- BESTmcmc(y, priors=priors, rnd.seed=123) )


