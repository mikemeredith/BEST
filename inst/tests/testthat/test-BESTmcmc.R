
# Tests for BESTmcmc and retro power with priors=NULL
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
  if(packageVersion("rjags") >= "4.0.0")  {
    expect_that(round(Bout2$mu1, 5),
      is_equivalent_to(c(4.98114, 4.88526, 4.25582, 4.79039, 4.71482, 4.70673, 4.80105, 4.75658, 5.14135)))
    expect_that(round(Bout2$mu2, 5),
      is_equivalent_to(c(3.19374, 3.28132, 3.17236, 3.56339, 3.49283, 3.44266, 3.10148, 3.25103, 3.36642)))
    expect_that(round(Bout2$nu, 5),
      is_equivalent_to(c(13.83349,  5.33221,  5.81619, 42.31951, 61.54992, 56.15759, 32.23504, 18.94366, 34.75114)))
    expect_that(round(Bout2$sigma1, 5),
      is_equivalent_to(c( 0.62254, 0.76634, 0.60601, 0.43419, 0.40904, 0.75910, 0.94848, 1.04851, 1.00506)))
    expect_that(round(Bout2$sigma2, 5),
      is_equivalent_to(c(1.35774, 2.27456, 1.85109, 0.47483, 0.73277, 0.47497, 0.51886, 0.55978, 0.89049)))
    expect_that(round(attr(Bout2, "Rhat"), 5),
      is_equivalent_to(c(1.06104,  2.67465,  4.24669,  2.91354,  3.57077)))
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
      is_equivalent_to(c(0.36364, 0.09091, 0.09091, 0.27273, 0.09091, 0.09091, 0.27273,
        0.27273, 0.18182, 0.09091, 0.09091, 0.09091)))
    expect_that(round(pow2[, 2], 5),
      is_equivalent_to(c(0.10678, 0.00000, 0.00000, 0.04645, 0.00000, 0.00000, 0.04645,
        0.04645, 0.00729, 0.00000, 0.00000, 0.00000)))
    expect_that(round(pow2[, 3], 5),
      is_equivalent_to(c(0.63320, 0.25887, 0.25887, 0.52243, 0.25887, 0.25887, 0.52243,
        0.52243, 0.39781, 0.25887, 0.25887, 0.25887)))
  }
})

y0 <- c(1.89, 1.78, 1.30, 1.74, 1.33, 0.89)
Bout1 <- BESTmcmc(y0, numSavedSteps = 9, burnInSteps = 1,
  verbose=FALSE, rnd.seed=123)

test_that("BESTmcmc with 1 group gives same output",  {
  expect_that(class(Bout1), equals(c("BEST", "data.frame")))
  expect_that(colnames(Bout1),
    equals(c("mu", "nu", "sigma")))
  if(packageVersion("rjags") >= "4.0.0")  {
    expect_that(round(Bout1$mu, 5),
      is_equivalent_to(c(1.31110, 1.42774, 1.31009, 1.49958, 1.66729, 1.66788, 1.52731, 1.51570, 1.44977)))
    expect_that(round(Bout1$nu, 5),
      is_equivalent_to(c(23.43798, 15.86208, 43.93298, 33.19201,  8.74280,  7.68395, 45.23169, 20.91182,  5.84763)))
    expect_that(round(Bout1$sigma, 5),
      is_equivalent_to(c(0.57750, 0.85844, 0.57488, 0.29597, 0.52261, 0.54746, 0.30884, 0.26675, 0.35829)))
    expect_that(round(attr(Bout1, "Rhat"), 5),
      is_equivalent_to(c(2.80791, 0.93063, 2.25087)))
    expect_that(attr(Bout1, "n.eff"),
      is_equivalent_to(c(9,  9,  9)))
  }
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
  if(packageVersion("rjags") >= "4.0.0")  {
    expect_that(round(pow1[, 1], 5),
      is_equivalent_to(c(0.81818, 0.09091, 0.09091, 0.72727, 0.09091, 0.09091,
        0.45455, 0.81818, 0.45455, 0.09091, 0.09091, 0.09091)))
    expect_that(round(pow1[, 2], 5),
      is_equivalent_to(c(0.60219, 0.00000, 0.00000, 0.47757, 0.00000, 0.00000,
        0.18155, 0.60219, 0.18155, 0.00000, 0.00000, 0.00000)))
    expect_that(round(pow1[, 3], 5),
      is_equivalent_to(c(0.99271, 0.25887, 0.25887, 0.95355, 0.25887, 0.25887,
        0.73161, 0.99271, 0.73161, 0.25887, 0.25887, 0.25887)))
  }
})
