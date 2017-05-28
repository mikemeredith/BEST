
# Tests for BESTmcmc

context("BESTmcmc_issues")

test_that("Old issues resolved with BESTmcmc",  {
  # Issue #6
  # silberzwiebel's data: should no throw error (but output is nonsense)
  y0 <- rep(-7:4, c(1, 2, 10, 12, 7, 31, 89, 231, 19, 4, 1, 1))
  expect_silent(BESTmcmc(y0, NULL, list(), numSavedSteps=100, burnInSteps = 100,
    parallel=FALSE, verbose=FALSE))

  # Issue #7 : R's global RNG set by BESTmcmc (rnd.seed should be passed directly to JAGS)
  # This was a jagsUI problem.
  set.seed(123) # I want reproducible samples
  means <- numeric(6)
  for(i in 1:6) {
    y0 <- rnorm(7)
    means[i] <- mean(y0)
    BESTmcmc(y0, NULL, list(), numSavedSteps=100, burnInSteps = 100,
      rnd.seed=(456), verbose=FALSE)
  }
  expect_false(any(diff(means) == 0)) # all samples should be different
})
