
# Tests for BESTmcmc

# library(BEST)
# library(testthat)
# test_file("test-BESTmcmc.R")

require(rjags)

context("BESTmcmc&retroPower")

y1 <- c(5.77, 5.33, 4.59, 4.33, 3.66, 4.48)
y2 <- c(3.88, 3.55, 3.29, 2.59, 2.33, 3.59)
Bout2 <- BESTmcmc(y1, y2, numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 2 groups gives same output",  {
  expect_that(class(Bout2), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout2),
    equals(c("mu1", "mu2", "nu", "sigma1", "sigma2")))
  expect_that(round(Bout2$mu1, 5), 
    is_equivalent_to(c(4.85621, 3.89109, 3.78930, 4.60468, 5.07111,
      4.10311, 4.61154, 4.80428, 4.30812)))
  expect_that(round(Bout2$mu2, 5), 
    is_equivalent_to(c(2.89325, 3.22934, 2.77639, 3.08105, 3.31025,
      3.03435, 2.97275, 2.87187, 3.01270)))
  expect_that(round(Bout2$nu, 5), 
    is_equivalent_to(c(14.34222,  24.07308,  21.08025,  70.13101,
      75.56236,  44.05460, 100.20045, 117.95162, 134.21188)))
  expect_that(round(Bout2$sigma1, 5), 
    is_equivalent_to(c(0.59383, 0.93246, 1.22044, 0.55669, 0.98399,
      1.25699, 0.83769, 0.62128, 0.96044)))
  expect_that(round(Bout2$sigma2, 5), 
    is_equivalent_to(c(1.07994, 0.81221, 0.91534, 1.31522, 0.57786,
      0.54844, 0.68362, 0.57276, 0.55713)))
  expect_that(round(attr(Bout2, "Rhat"), 5), 
    is_equivalent_to(c(1.10123,  1.21913,  5.21131,  0.89490,  1.22130)))
  expect_that(attr(Bout2, "n.eff"), 
    is_equivalent_to(c(9,  9,  9,  9,  9)))
})

test_that("BESTpower retro with 2 groups gives same output",  {
  pow2 <- BESTpower(Bout2, 
    ROPEm=c(-0.1,0.1), ROPEsd=c(-2,2), ROPEeff=c(-0.5,0.5), 
    maxHDIWm=2.0, maxHDIWsd=2.0, maxHDIWeff=2.0,
    nRep=9, mcmcLength=1000, verbose=FALSE, rnd.seed=456)
  expect_that(class(pow2), equals("matrix"))
  expect_that(colnames(pow2),
    equals(c("mean", "CrIlo", "CrIhi")))
  expect_that(rownames(pow2),
    equals(c("  mean:   HDI > ROPE", "  mean:   HDI < ROPE",
      "  mean:  HDI in ROPE", "  mean: HDI width ok",
      "    sd:   HDI > ROPE", "    sd:   HDI < ROPE",
      "    sd:  HDI in ROPE", "    sd: HDI width ok",
      "effect:   HDI > ROPE", "effect:   HDI < ROPE",
      "effect:  HDI in ROPE", "effect: HDI width ok")))
  expect_that(round(pow2[, 1], 5), 
    is_equivalent_to(c(0.54545, 0.09091, 0.09091, 0.18182, 0.09091,
      0.09091, 0.63636, 0.27273, 0.36364, 0.09091, 0.09091, 0.09091)))
  expect_that(round(pow2[, 2], 5), 
    is_equivalent_to(c(0.26839, 0.00000, 0.00000, 0.00729, 0.00000,
      0.00000, 0.36680, 0.04645, 0.10678, 0.00000, 0.00000, 0.00000)))
  expect_that(round(pow2[, 3], 5), 
    is_equivalent_to(c(0.81845, 0.25887, 0.25887, 0.39781, 0.25887,
      0.25887, 0.89322, 0.52243, 0.63320, 0.25887, 0.25887, 0.25887)))
})

y0 <- c(1.89, 1.78, 1.30, 1.74, 1.33, 0.89)
Bout1 <- BESTmcmc(y0, numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 1 group gives same output",  {
  expect_that(class(Bout1), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout1),
    equals(c("mu", "nu", "sigma")))
  expect_that(round(Bout1$mu, 5), 
    is_equivalent_to(c(1.58291, 1.57144, 1.57846, 1.51096, 1.33814,
      1.34559, 1.40930, 1.49599, 1.51241)))
  expect_that(round(Bout1$nu, 5), 
    is_equivalent_to(c(30.24444, 42.84582, 41.61960, 35.27966,
      34.31609, 10.08930, 22.43569, 26.38340, 20.31474)))
  expect_that(round(Bout1$sigma, 5), 
    is_equivalent_to(c(0.36703, 0.34621, 0.26857, 0.48060, 0.74017,
      0.92785, 0.40013, 0.33565, 0.25275)))
  expect_that(round(attr(Bout1, "Rhat"), 5), 
    is_equivalent_to(c(2.19598, 1.45654, 2.53052)))
  expect_that(attr(Bout1, "n.eff"), 
    is_equivalent_to(c(9,  9,  9)))
})

test_that("BESTpower retro with 1 group gives same output",  {
  pow1 <- BESTpower(Bout1, 
    ROPEm=c(-0.5,0.5), ROPEsd=c(-1,1), ROPEeff=c(-1,1), 
    maxHDIWm=2.0, maxHDIWsd=2.0, maxHDIWeff=2.0,
    nRep=9, mcmcLength=1000, verbose=FALSE, rnd.seed=456)
  expect_that(class(pow1), equals("matrix"))
  expect_that(colnames(pow1),
    equals(c("mean", "CrIlo", "CrIhi")))
  expect_that(rownames(pow1),
    equals(c("  mean:   HDI > ROPE", "  mean:   HDI < ROPE",
      "  mean:  HDI in ROPE", "  mean: HDI width ok",
      "    sd:   HDI > ROPE", "    sd:   HDI < ROPE",
      "    sd:  HDI in ROPE", "    sd: HDI width ok",
      "effect:   HDI > ROPE", "effect:   HDI < ROPE",
      "effect:  HDI in ROPE", "effect: HDI width ok")))
  expect_that(round(pow1[, 1], 5), 
    is_equivalent_to(c(0.72727, 0.09091, 0.09091, 0.72727, 0.09091,
      0.09091, 0.54545, 0.81818, 0.54545, 0.09091, 0.09091, 0.09091)))
  expect_that(round(pow1[, 2], 5), 
    is_equivalent_to(c(0.47757, 0.00000, 0.00000, 0.47757, 0.00000,
      0.00000, 0.26839, 0.60219, 0.26839, 0.00000, 0.00000, 0.00000)))
  expect_that(round(pow1[, 3], 5), 
    is_equivalent_to(c(0.95355, 0.25887, 0.25887, 0.95355, 0.25887,
      0.25887, 0.81845, 0.99271, 0.81845, 0.25887, 0.25887, 0.25887)))
})
  