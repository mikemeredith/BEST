
# Tests for BESTmcmc and retro power with gamma priors
#  (both work with the same BESTmcmc output object).


context("BESTmcmc&retroPower_gammaPriors")

y1 <- c(5.77, 5.33, 4.59, 4.33, 3.66, 4.48)
y2 <- c(3.88, 3.55, 3.29, 2.59, 2.33, 3.59)
Bout2 <- BESTmcmc(y1, y2, priors=list(),
  numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 2 groups and default gamma priors gives same output",  {
  expect_that(class(Bout2), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout2),
    equals(c("mu1", "mu2", "nu", "sigma1", "sigma2")))
  if(packageVersion("rjags") >= "4.0.0")  {
    expect_that(round(Bout2$mu1, 5),
      is_equivalent_to(c(4.61437, 4.44746, 4.35208, 5.19768, 5.24549, 4.85071, 4.48226, 4.78731, 3.19296)))
    expect_that(round(Bout2$mu2, 5),
      is_equivalent_to(c(3.38234, 3.31267, 3.25178, 3.44527, 3.42282, 3.29180, 3.21568, 3.44653, 2.98123)))
    expect_that(round(Bout2$nu, 5),
      is_equivalent_to(c(21.94517, 13.35176,  7.70370, 31.38732,  5.16274, 26.68239, 33.93732, 27.87325, 30.40717)))
    expect_that(round(Bout2$sigma1, 5),
      is_equivalent_to(c(1.03883, 1.01869, 1.12193, 0.81696, 0.85216, 0.86884, 0.70153, 1.40755, 1.96448)))
    expect_that(round(Bout2$sigma2, 5),
      is_equivalent_to(c(1.02930, 0.97967, 1.07990, 0.53815, 0.54342, 0.63613, 0.66597, 0.66632, 0.58443)))
    expect_that(round(attr(Bout2, "Rhat"), 5),
      is_equivalent_to(c(1.62124,  1.23876,  1.49260,  1.38925,  7.26812 )))
    expect_that(attr(Bout2, "n.eff"),
      is_equivalent_to(c(9,  9,  9,  9,  9)))
  }
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
  if(packageVersion("rjags") >= "4.0.0")  {
    expect_that(round(pow2[, 1], 5),
      is_equivalent_to(c(0.27273, 0.09091, 0.09091, 0.18182, 0.09091, 0.09091, 
        0.36364, 0.18182, 0.18182, 0.09091, 0.09091, 0.09091)))
    expect_that(round(pow2[, 2], 5),
      is_equivalent_to(c(0.04645, 0.00000, 0.00000, 0.00729, 0.00000, 0.00000,
        0.10678, 0.00729, 0.00729, 0.00000, 0.00000, 0.00000)))
    expect_that(round(pow2[, 3], 5),
      is_equivalent_to(c(0.52243, 0.25887, 0.25887, 0.39781, 0.25887, 0.25887,
        0.63320, 0.39781, 0.39781, 0.25887, 0.25887, 0.25887)))
  }
})

Bout2a <- BESTmcmc(y1, y2,
  priors=list(muM=7:8, muSD=10:11, sigmaMode=4, sigmaSD=8, nuMean=5, nuSD=10),
  numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 2 groups and informative gamma priors gives same output",  {
  expect_that(class(Bout2a), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout2a),
    equals(c("mu1", "mu2", "nu", "sigma1", "sigma2")))
if(packageVersion("rjags") >= "4.0.0")  {
   expect_that(round(Bout2a$mu1, 5),
      is_equivalent_to(c(4.93248, 4.86504, 4.55111, 5.12664, 4.60494, 4.53156, 5.20077, 5.07907, 3.94995)))
    expect_that(round(Bout2a$mu2, 5),
      is_equivalent_to(c(3.08729, 3.14780, 3.41774, 3.98207, 3.87674, 3.64305, 3.25928, 3.20982, 3.09678)))
    expect_that(round(Bout2a$nu, 5),
      is_equivalent_to(c(31.18633, 33.74027, 26.64342, 26.20644, 38.20690, 43.69588,  5.91735,  9.26451, 19.25745)))
    expect_that(round(Bout2a$sigma1, 5),
      is_equivalent_to(c(0.63443, 1.03655, 0.77387, 1.28654, 1.10640, 1.70681, 0.54008, 1.04993, 1.57496)))
    expect_that(round(Bout2a$sigma2, 5),
      is_equivalent_to(c(0.90291, 0.75760, 0.85716, 1.00634, 0.85724, 0.96690, 0.49994, 0.65312, 0.52598)))
    expect_that(round(attr(Bout2a, "Rhat"), 5),
      is_equivalent_to(c(0.95010,  3.68356,  2.89774,  1.32351,  3.86606)))
    expect_that(attr(Bout2a, "n.eff"),
      is_equivalent_to(c(9,  9,  9,  9,  9)))
  }
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
if(packageVersion("rjags") >= "4.0.0")  {
    expect_that(round(Bout1$mu, 5),
      is_equivalent_to(c(1.55073, 1.37981, 1.37126, 1.54232, 1.42041, 1.64843, 1.42762, 1.53105, 1.39775)))
    expect_that(round(Bout1$nu, 5),
      is_equivalent_to(c(32.06687, 89.00384, 86.96579, 20.82893,  5.57220,  3.98534, 45.69073, 43.44243, 30.36267)))
    expect_that(round(Bout1$sigma, 5),
      is_equivalent_to(c(0.43361, 0.58152, 0.57456, 0.31313, 0.27469, 0.45573, 0.53270, 0.26520, 0.36109)))
    expect_that(round(attr(Bout1, "Rhat"), 5),
      is_equivalent_to(c(1.14463, 2.34712, 1.50120 )))
    expect_that(attr(Bout1, "n.eff"),
      is_equivalent_to(c(9,  9,  9)))
  }
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
if(packageVersion("rjags") >= "4.0.0")  {
    expect_that(round(pow1[, 1], 5),
      is_equivalent_to(c( 0.90909, 0.09091, 0.09091, 0.90909, 0.09091, 0.09091,
        0.45455, 0.90909, 0.45455, 0.09091, 0.09091, 0.09091)))
    expect_that(round(pow1[, 2], 5),
      is_equivalent_to(c(0.74113, 0.00000, 0.00000, 0.74113, 0.00000, 0.00000,
        0.18155, 0.74113, 0.18155, 0.00000, 0.00000, 0.00000)))
    expect_that(round(pow1[, 3], 5),
      is_equivalent_to(c(1.00000, 0.25887, 0.25887, 1.00000, 0.25887, 0.25887, 
        0.73161, 1.00000, 0.73161, 0.25887, 0.25887, 0.25887)))
  }
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
if(packageVersion("rjags") >= "4.0.0")  {
    expect_that(round(Bout1a$mu, 5),
      is_equivalent_to(c(1.68589, 1.66578, 1.75534, 1.56379, 1.62999, 1.28949, 1.61038, 1.38212, 1.60525)))
    expect_that(round(Bout1a$nu, 5),
      is_equivalent_to(c(1.82353, 6.20606, 5.21446, 1.75047, 2.27062, 2.49363, 3.09216, 7.59883, 3.83828)))
    expect_that(round(Bout1a$sigma, 5),
      is_equivalent_to(c( 0.37997, 0.36664, 0.49527, 0.61515, 0.58572, 0.67730, 0.41813, 0.35117, 0.42800)))
    expect_that(round(attr(Bout1a, "Rhat"), 5),
      is_equivalent_to(c(1.44123, 1.30808, 3.53620 )))
    expect_that(attr(Bout1a, "n.eff"),
      is_equivalent_to(c(9,  9,  9)))
  }
  expect_that(names(attr(Bout1a, "priors")), equals(c("muM", "muSD", "sigmaMode",
    "sigmaSD", "nuMean", "nuSD")))
})
