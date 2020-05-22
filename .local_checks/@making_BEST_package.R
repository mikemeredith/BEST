# Individual checks

# setwd("D:/Github/BEST_package")
setwd("../..")

library(devtools)
sIg <- scan("spellcheckIgnore.txt", what='character', comment="#")
tmp <- spell_check("BEST", ignore=c(hunspell::en_stats, sIg), "en_GB")
length(tmp)  # number of misspellings found
tmp  # error if length == 0

devtools::load_all("C:/GitHub/BEST_package/BEST")

unlink(list.files(pattern="Rplots.pdf", recursive=TRUE))
system("R CMD build BEST")  # Produces the .tar.gz file
system("R CMD check BEST_0.5.2.9000.tar.gz")
system("R CMD check --run-donttest BEST_0.5.2.9000.tar.gz")
# system("R CMD check --as-cran BEST_0.5.2.9000.tar.gz --no-manual")
system("R CMD INSTALL --build BEST_0.5.2.9000.tar.gz") # installs and produces the .zip binary
# system("R CMD INSTALL BEST_0.5.2.9000.tar.gz") # installs only

system("R CMD INSTALL BEST") # Use this for a "dev" install.

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


