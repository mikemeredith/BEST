
# Tests for BESTmcmc and retro power with gamma priors
#  (both work with the same BESTmcmc output object).


context("BESTmcmc&retroPower_gammaPriors")

y1 <- c(5.77, 5.33, 4.59, 4.33, 3.66, 4.48)
y2 <- c(3.88, 3.55, 3.29, 2.59, 2.33, 3.59)
Bout2s <- BESTmcmc(y1, y2, priors=list(),
  numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123, parallel=FALSE)
Bout2p <- BESTmcmc(y1, y2, priors=list(),
  numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123, parallel=TRUE)

test_that("BESTmcmc with 2 groups and default gamma priors gives same output",  {
  expect_equivalent(Bout2s, Bout2p)
  expect_that(class(Bout2s), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout2s),
    equals(c("mu1", "mu2", "nu", "sigma1", "sigma2")))
  if(packageVersion("rjags") >= "4.0.0")  {
    expect_equivalent(round(colMeans(Bout2s), 5),
        c(4.50176,  3.00586, 38.12665,  0.78138,  0.99560))
  }
})

test_that("BESTpower retro with 2 groups gives same output",  {
  pow2s <- BESTpower(Bout2s,
    ROPEm=c(-0.1,0.1), ROPEsd=c(-2,2), ROPEeff=c(-0.5,0.5),
    maxHDIWm=2.0, maxHDIWsd=2.0, maxHDIWeff=2.0,
    nRep=9, mcmcLength=1000, verbose=FALSE, rnd.seed=456, parallel=FALSE)
  pow2p <- BESTpower(Bout2s,
    ROPEm=c(-0.1,0.1), ROPEsd=c(-2,2), ROPEeff=c(-0.5,0.5),
    maxHDIWm=2.0, maxHDIWsd=2.0, maxHDIWeff=2.0,
    nRep=9, mcmcLength=1000, verbose=FALSE, rnd.seed=456, parallel=TRUE)
  expect_equivalent(pow2s, pow2p)
  expect_that(class(pow2s), equals("matrix"))
  expect_that(colnames(pow2s),
    equals(c("mean", "CrIlo", "CrIhi")))
  expect_that(rownames(pow2s),
    equals(c("  mean:   HDI > ROPE", "  mean:   HDI < ROPE",
      "  mean:  HDI in ROPE", "  mean: HDI width ok",
      "    sd:   HDI > ROPE", "    sd:   HDI < ROPE",
      "    sd:  HDI in ROPE", "    sd: HDI width ok",
      "effect:   HDI > ROPE", "effect:   HDI < ROPE",
      "effect:  HDI in ROPE", "effect: HDI width ok")))
  if(packageVersion("rjags") >= "4.0.0")  {
    expect_equivalent(round(colMeans(pow2s), 5), c(0.18182, 0.03635, 0.38181))
  }
})

Bout2a <- BESTmcmc(y1, y2,
  priors=list(muM=7:8, muSD=10:11, sigmaMode=4, sigmaSD=8, nuMean=5, nuSD=10),
  numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 2 groups and informative gamma priors gives same output",  {
  expect_equal(class(Bout2a), c("BEST", "data.frame"))
  expect_that(colnames(Bout2a),
    equals(c("mu1", "mu2", "nu", "sigma1", "sigma2")))
if(packageVersion("rjags") >= "4.0.0")  {
   expect_equivalent(round(colMeans(Bout2a), 5),
      c(4.70158, 3.47429, 7.12125, 0.99167, 0.71483))
    expect_equivalent(round(mean(attr(Bout2a, "Rhat")), 5), 1.83177)
    expect_equivalent(attr(Bout2a, "n.eff"), c(9,  9,  9,  9,  9))
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
    expect_equivalent(round(colMeans(Bout1), 5), c(1.53659, 21.29717,  0.36009))
    expect_equivalent(round(mean(attr(Bout1, "Rhat")), 5), 2.18427)
    expect_equivalent(attr(Bout1, "n.eff"), c(9,  9,  9))
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
    expect_equivalent(round(colMeans(pow1), 5), c(0.34848, 0.20931, 0.51474))
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
    expect_equivalent(round(colMeans(Bout1a), 5), c(1.53828, 5.78369, 0.38535))
    expect_equivalent(round(mean(attr(Bout1a, "Rhat")), 5), 1.33725)
    expect_that(attr(Bout1a, "n.eff"),
      is_equivalent_to(c(9,  9,  9)))
  }
  expect_that(names(attr(Bout1a, "priors")), equals(c("muM", "muSD", "sigmaMode",
    "sigmaSD", "nuMean", "nuSD")))
})
