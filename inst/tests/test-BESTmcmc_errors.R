
# Tests for BESTmcmc 

context("BESTmcmc_errors")

y1 <- c(5.77, 5.33, 4.59, 4.33, 3.66, 4.48)
y2 <- c(3.88, 3.55, 3.29, 2.59, 2.33, 3.59)
Bout2 <- BESTmcmc(y1, y2, numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc gives sensible errors",  {
  expect_error(BESTmcmc(1:2, 1:4, numSavedSteps = 9, burnInSteps = 1),
    "Minimum size for both samples is 3")
  expect_error(BESTmcmc(c(1:4, NA), 1:4, numSavedSteps = 9, burnInSteps = 1),
    "The input data include NA or Inf")
  expect_error(BESTmcmc(1:4, 1:4, 9, 1),
    "'priors' is now the 3rd argument")
  expect_error(BESTmcmc(1:4, 1:4, 9, 1),
    "it must be a list")
  expect_error(BESTmcmc(1:4, 1:4, priors="A", 9, 1),
    "'priors' must be a list")
  expect_error(BESTmcmc(1:4, 1:4, list(nonsense=8), 9, 1),
    "Invalid items in prior specification")
  expect_error(BESTmcmc(1:4, 1:4, list(muSD=-1), 9, 1),
    "muSD must be > 0")
  expect_error(BESTmcmc(1:4, 1:4, list(sigmaMode=-1), 9, 1),
    "gamma prior must be > 0")
  expect_error(BESTmcmc(1:4, 1:4, list(nuSD=-1), 9, 1),
    "gamma prior must be > 0")

})