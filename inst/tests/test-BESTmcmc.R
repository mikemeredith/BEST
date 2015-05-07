
# Tests for BESTmcmc and retro power
#  (both work with the same BESTmcmc output object).


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
    is_equivalent_to(c(4.58957, 4.63928, 4.56552, 4.68522, 4.79811,
      4.87478, 4.31500, 4.88270, 4.11194)))
  expect_that(round(Bout2$mu2, 5), 
    is_equivalent_to(c(3.30069, 3.32236, 3.12548, 3.27334, 3.50424,
      2.98577, 3.05609, 3.44276, 3.16167)))
  expect_that(round(Bout2$nu, 5), 
    is_equivalent_to(c(5.35834, 7.88822, 12.93882, 32.09664, 31.80601,
      17.1242, 29.13672, 36.2641, 36.17713)))
  expect_that(round(Bout2$sigma1, 5), 
    is_equivalent_to(c(0.38631, 1.0764, 0.98101, 0.52098, 0.80038, 0.73821,
      0.69582, 0.91066, 0.94981)))
  expect_that(round(Bout2$sigma2, 5), 
    is_equivalent_to(c(0.48743, 1.07449, 1.23421, 0.51216, 0.56319, 0.48284,
      0.95822, 0.93705, 0.77755)))
  expect_that(round(attr(Bout2, "Rhat"), 5), 
    is_equivalent_to(c(1.40241, 0.88268, 3.33835, 1.01467, 1.62462)))
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
    is_equivalent_to(c(0.45455, 0.09091, 0.09091, 0.36364, 0.09091, 0.09091,
      0.54545, 0.36364, 0.27273, 0.09091, 0.09091, 0.09091)))
  expect_that(round(pow2[, 2], 5), 
    is_equivalent_to(c(0.18155, 0, 0, 0.10678, 0, 0, 0.26839,
      0.10678, 0.04645, 0, 0, 0)))
  expect_that(round(pow2[, 3], 5), 
    is_equivalent_to(c(0.73161, 0.25887, 0.25887, 0.6332, 0.25887,
      0.25887, 0.81845, 0.6332, 0.52243, 0.25887, 0.25887, 0.25887)))
})

y0 <- c(1.89, 1.78, 1.30, 1.74, 1.33, 0.89)
Bout1 <- BESTmcmc(y0, numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 1 group gives same output",  {
  expect_that(class(Bout1), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout1),
    equals(c("mu", "nu", "sigma")))
  expect_that(round(Bout1$mu, 5), 
    is_equivalent_to(c(1.53373, 1.43131, 1.72584, 1.98664, 1.96624,
      1.95121, 1.49141, 1.68495, 1.58938)))
  expect_that(round(Bout1$nu, 5), 
    is_equivalent_to(c(37.64804, 45.35948, 1.39458, 10.80727, 5.86184,
      26.89274, 5.43921, 13.67402, 11.39325)))
  expect_that(round(Bout1$sigma, 5), 
    is_equivalent_to(c(0.43763, 0.4345, 0.34226, 0.58617, 0.6972, 0.67925,
      0.56161, 0.48649, 0.55028)))
  expect_that(round(attr(Bout1, "Rhat"), 5), 
    is_equivalent_to(c(3.33436, 1.2654, 3.64821)))
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
    is_equivalent_to(c(0.81818, 0.09091, 0.09091, 0.90909, 0.09091,
      0.09091, 0.54545, 0.90909, 0.45455, 0.09091, 0.09091, 0.09091)))
  expect_that(round(pow1[, 2], 5), 
    is_equivalent_to(c(0.60219, 0, 0, 0.74113, 0, 0, 0.26839, 0.74113,
      0.18155, 0, 0, 0)))
  expect_that(round(pow1[, 3], 5), 
    is_equivalent_to(c(0.99271, 0.25887, 0.25887, 1, 0.25887, 0.25887,
      0.81845, 1, 0.73161, 0.25887, 0.25887, 0.25887)))
})
  