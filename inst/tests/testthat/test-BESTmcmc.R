
# Tests for BESTmcmc and retro power with priors=NULL
#  (both work with the same BESTmcmc output object).


context("BESTmcmc&retroPower")

y1 <- c(5.77, 5.33, 4.59, 4.33, 3.66, 4.48)
y2 <- c(3.88, 3.55, 3.29, 2.59, 2.33, 3.59)
Bout2s <- BESTmcmc(y1, y2, numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123, parallel=FALSE)
Bout2p <- BESTmcmc(y1, y2, numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123, parallel=TRUE)

test_that("BESTmcmc with 2 groups gives same output",  {
  expect_that(class(Bout2s), equals(c("BEST", "data.frame")))
  expect_that(class(Bout2p), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout2s),
    equals(c("mu1", "mu2", "nu", "sigma1", "sigma2")))
  expect_that(colnames(Bout2p),
    equals(c("mu1", "mu2", "nu", "sigma1", "sigma2")))
  expect_equivalent(Bout2s, Bout2p)
  if(packageVersion("rjags") >= "4.0.0")  {
    expect_equivalent(round(colMeans(Bout2s), 5),
      c(4.73030,  3.06054, 28.10154,  0.77949,  0.69464))
    expect_equal(round(mean(attr(Bout2s, "Rhat")), 5), 1.28783)
    expect_equal(round(mean(attr(Bout2p, "Rhat")), 5), round(mean(attr(Bout2s, "Rhat")), 5))
    expect_equivalent(attr(Bout2s, "n.eff"), c(9,  9,  9,  9,  9))
    expect_equivalent(attr(Bout2p, "n.eff"), c(9,  9,  9,  9,  9))
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
  expect_equal(pow2s, pow2p)
  expect_that(class(pow2s), equals(c("matrix", "array")))
  expect_that(colnames(pow2s),
    equals(c("mean", "CrIlo", "CrIhi")))
  expect_that(rownames(pow2s),
    equals(c("  mean:   HDI > ROPE", "  mean:   HDI < ROPE",
      "  mean:  HDI in ROPE", "  mean: HDI width ok",
      "    sd:   HDI > ROPE", "    sd:   HDI < ROPE",
      "    sd:  HDI in ROPE", "    sd: HDI width ok",
      "effect:   HDI > ROPE", "effect:   HDI < ROPE",
      "effect:  HDI in ROPE", "effect: HDI width ok")))
  expect_equivalent(round(colMeans(pow2s), 5), c(0.26515, 0.12176, 0.45393))
})

y0 <- c(1.89, 1.78, 1.30, 1.74, 1.33, 0.89)
Bout1s <- BESTmcmc(y0, numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123, parallel=FALSE)
Bout1p <- BESTmcmc(y0, numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123, parallel=TRUE)

test_that("BESTmcmc with 1 group gives same output",  {
  expect_that(class(Bout1s), equals(c("BEST", "data.frame")))
  expect_that(class(Bout1p), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout1s),
    equals(c("mu", "nu", "sigma")))
  expect_that(colnames(Bout1s),
    equals(c("mu", "nu", "sigma")))
  expect_equivalent(round(colMeans(Bout1s), 5), c(1.39340, 24.30192,  0.48864))
  expect_equal(round(mean(attr(Bout1s, "Rhat")), 5), 1.36887)
  expect_equal(round(mean(attr(Bout1p, "Rhat")), 5), round(mean(attr(Bout1s, "Rhat")), 5))
  expect_equal(attr(Bout1p, "Rhat"), attr(Bout1s, "Rhat"))
  expect_equivalent(attr(Bout1s, "n.eff"), c(9,  9,  9))
  expect_equivalent(attr(Bout1p, "n.eff"), c(9,  9,  9))
})

test_that("BESTpower retro with 1 group gives same output",  {
  pow1s <- BESTpower(Bout1s,
    ROPEm=c(-0.5,0.5), ROPEsd=c(-1,1), ROPEeff=c(-1,1),
    maxHDIWm=2.0, maxHDIWsd=2.0, maxHDIWeff=2.0,
    nRep=9, mcmcLength=1000, verbose=FALSE, rnd.seed=456, parallel=FALSE)
  pow1p <- BESTpower(Bout1s,
    ROPEm=c(-0.5,0.5), ROPEsd=c(-1,1), ROPEeff=c(-1,1),
    maxHDIWm=2.0, maxHDIWsd=2.0, maxHDIWeff=2.0,
    nRep=9, mcmcLength=1000, verbose=FALSE, rnd.seed=456, parallel=TRUE)
  expect_equal(pow1s, pow1p)
  expect_equal(class(pow1s), c("matrix", "array"))
  expect_that(colnames(pow1s),
    equals(c("mean", "CrIlo", "CrIhi")))
  expect_that(rownames(pow1s),
    equals(c("  mean:   HDI > ROPE", "  mean:   HDI < ROPE",
      "  mean:  HDI in ROPE", "  mean: HDI width ok",
      "    sd:   HDI > ROPE", "    sd:   HDI < ROPE",
      "    sd:  HDI in ROPE", "    sd: HDI width ok",
      "effect:   HDI > ROPE", "effect:   HDI < ROPE",
      "effect:  HDI in ROPE", "effect: HDI width ok")))
  expect_equivalent(round(colMeans(pow1s), 5), c(0.28030, 0.13426, 0.46934 ))
})
