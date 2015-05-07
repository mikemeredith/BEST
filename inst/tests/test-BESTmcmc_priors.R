
# Tests for BESTmcmc and retro power with informative priors
#  (both work with the same BESTmcmc output object).


context("BESTmcmc&retroPower_priors")

y1 <- c(5.77, 5.33, 4.59, 4.33, 3.66, 4.48)
y2 <- c(3.88, 3.55, 3.29, 2.59, 2.33, 3.59)
Bout2 <- BESTmcmc(y1, y2, priors=list(),
  numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 2 groups and default gamma priors gives same output",  {
  expect_that(class(Bout2), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout2),
    equals(c("mu1", "mu2", "nu", "sigma1", "sigma2")))
  expect_that(round(Bout2$mu1, 5),
    is_equivalent_to(c(4.53926, 4.61859, 4.57063, 4.29667, 3.99112,
      4.66591, 4.76619, 4.48431, 4.38613)))
  expect_that(round(Bout2$mu2, 5),
    is_equivalent_to(c(3.2797, 3.54284, 3.0141, 2.9892, 2.99361,
      2.50182, 3.4123, 3.4485, 3.26311)))
  expect_that(round(Bout2$nu, 5),
    is_equivalent_to(c(23.82615, 26.2636, 5.64102, 48.69661, 55.59175,
      42.32349, 67.19148, 70.25186, 161.87068)))
  expect_that(round(Bout2$sigma1, 5),
    is_equivalent_to(c(0.74105, 0.68235, 0.88782, 1.55269, 1.52947,
      1.39599, 0.65826, 0.69805, 0.77205)))
  expect_that(round(Bout2$sigma2, 5),
    is_equivalent_to(c(0.66695, 0.79608, 0.84857, 0.70704, 1.4542,
      0.90285, 0.5798, 0.58487, 0.69581)))
  expect_that(round(attr(Bout2, "Rhat"), 5),
    is_equivalent_to(c(1.25864, 2.02649, 2.09094, 7.67192, 1.5592)))
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
    is_equivalent_to(c(0.54545, 0.09091, 0.09091, 0.27273, 0.09091,
      0.09091, 0.45455, 0.36364, 0.18182, 0.09091, 0.09091, 0.09091)))
  expect_that(round(pow2[, 2], 5),
    is_equivalent_to(c(0.26839, 0, 0, 0.04645, 0, 0, 0.18155, 0.10678,
      0.00729, 0, 0, 0)))
  expect_that(round(pow2[, 3], 5),
    is_equivalent_to(c(0.81845, 0.25887, 0.25887, 0.52243, 0.25887,
      0.25887, 0.73161, 0.63320, 0.39781, 0.25887, 0.25887, 0.25887)))
})

Bout2a <- BESTmcmc(y1, y2,
  priors=list(muM=7:8, muSD=10:11, sigmaMode=4, sigmaSD=8, nuMean=5, nuSD=10),
  numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 2 groups and informative gamma priors gives same output",  {
  expect_that(class(Bout2a), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout2a),
    equals(c("mu1", "mu2", "nu", "sigma1", "sigma2")))
  expect_that(round(Bout2a$mu1, 5),
    is_equivalent_to(c(4.50942, 4.76061, 4.55091, 4.23131, 5.20829,
      4.67565, 5.17247, 5.33514, 5.0565)))
  expect_that(round(Bout2a$mu2, 5),
    is_equivalent_to(c(3.12087, 3.59951, 3.51425, 3.29169, 3.3975,
      3.43153, 3.19523, 3.13906, 3.37815)))
  expect_that(round(Bout2a$nu, 5),
    is_equivalent_to(c(6.70917, 4.0176, 6.0169, 3.24085, 9.41055,
      7.48556, 13.39461, 8.19616, 14.11822)))
  expect_that(round(Bout2a$sigma1, 5),
    is_equivalent_to(c(0.74704, 1.2973, 1.39195, 1.24604, 0.96313,
      0.6738, 0.83346, 1.65861, 0.86389)))
  expect_that(round(Bout2a$sigma2, 5),
    is_equivalent_to(c(0.36881, 1.06708, 0.97023, 0.60223, 0.61518, 0.9242, 0.66605, 0.77598, 0.48513)))
  expect_that(round(attr(Bout2a, "Rhat"), 5),
    is_equivalent_to(c(1.67509,  1.15600,  1.98096,  0.89778, 1.00064)))
  expect_that(attr(Bout2a, "n.eff"),
    is_equivalent_to(c(9,  9,  9,  9,  9)))
  PR <- attr(Bout2a, "priors")
  expect_that(names(PR), equals(c("muM", "muSD", "sigmaMode", "sigmaSD",
    "nuMean", "nuSD")))
  expect_that(PR$muM, equals(7:8))
  expect_that(PR$muSD, equals(10:11))
})


#### One group tests
#### ===============
y0 <- c(1.89, 1.78, 1.30, 1.74, 1.33, 0.89)
Bout1 <- BESTmcmc(y0, priors=list(),
  numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 1 group and default gamma priors gives same output",  {
  expect_that(class(Bout1), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout1),
    equals(c("mu", "nu", "sigma")))
  expect_that(round(Bout1$mu, 5),
    is_equivalent_to(c(1.55405, 1.63076, 1.83076, 1.54219, 1.46405,
      1.5513, 1.54813, 1.93568, 2.09194)))
  expect_that(round(Bout1$nu, 5),
    is_equivalent_to(c(12.45647, 19.08607, 22.53563, 95.73921, 9.64915,
      13.33534, 11.49823, 14.34836, 36.91752)))
  expect_that(round(Bout1$sigma, 5),
    is_equivalent_to(c(0.68282, 0.5537, 0.72221, 0.27075, 0.30985,
      0.26813, 0.56542, 0.57178, 0.61748)))
  expect_that(round(attr(Bout1, "Rhat"), 5),
    is_equivalent_to(c(1.59326, 1.13276, 5.36971)))
  expect_that(attr(Bout1, "n.eff"),
    is_equivalent_to(c(9,  9,  9)))
})

test_that("BESTpower retro with 1 group and default gamma priors gives same output",  {
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
    is_equivalent_to(c(0.81818, 0.09091, 0.09091, 0.81818, 0.09091,
      0.09091, 0.45455, 0.81818, 0.54545, 0.09091, 0.09091, 0.09091)))
  expect_that(round(pow1[, 2], 5),
    is_equivalent_to(c(0.60219, 0, 0, 0.60219, 0, 0, 0.18155, 0.60219,
      0.26839, 0, 0, 0)))
  expect_that(round(pow1[, 3], 5),
    is_equivalent_to(c(0.99271, 0.25887, 0.25887, 0.99271, 0.25887,
      0.25887, 0.73161, 0.99271, 0.81845, 0.25887, 0.25887, 0.25887)))
})


y0 <- c(1.89, 1.78, 1.30, 1.74, 1.33, 0.89)
Bout1a <- BESTmcmc(y0,
  priors=list(muM=2, muSD=10, sigmaMode=3, sigmaSD=12, nuMean=5, nuSD=10),
  numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 1 group and informative gamma priors gives same output",  {
  expect_that(class(Bout1a), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout1a),
    equals(c("mu", "nu", "sigma")))
  expect_that(round(Bout1a$mu, 5),
    is_equivalent_to(c(1.37523, 1.01679, 0.92417, 1.23077, 1.54277,
      1.56209, 1.57113, 1.55808, 1.63596)))
  expect_that(round(Bout1a$nu, 5),
    is_equivalent_to(c(2.80934, 1.28834, 1.55243, 15.25831, 14.81094,
      4.64692, 21.73779, 16.51111, 7.03579)))
  expect_that(round(Bout1a$sigma, 5),
    is_equivalent_to(c(0.55073, 0.53202, 0.75882, 0.36801, 0.91949,
      0.40246, 0.39679, 0.31751, 0.22631)))
  expect_that(round(attr(Bout1a, "Rhat"), 5),
    is_equivalent_to(c(2.22914, 1.96134, 1.40245 )))
  expect_that(attr(Bout1a, "n.eff"),
    is_equivalent_to(c(9,  9,  9)))
  expect_that(names(attr(Bout1a, "priors")), equals(c("muM", "muSD", "sigmaMode",
    "sigmaSD", "nuMean", "nuSD")))
})
