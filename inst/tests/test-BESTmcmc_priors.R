
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
    is_equivalent_to(c(4.42717, 4.45490, 4.18039, 4.94525, 4.74881,
      4.49847, 4.01811, 3.69698, 3.87005)))
  expect_that(round(Bout2$mu2, 5),
    is_equivalent_to(c(2.47140, 3.55410, 3.09972, 3.05757, 3.09154,
      2.70293, 2.28640, 2.30421, 2.93424)))
  expect_that(round(Bout2$nu, 5),
    is_equivalent_to(c(30.75659, 20.20782, 49.19266, 17.11728, 18.83133,
      8.76044, 12.81865,  9.91033, 8.21868)))
  expect_that(round(Bout2$sigma1, 5),
    is_equivalent_to(c(1.55756, 0.89037, 0.91713, 0.58793, 0.49769, 0.97293,
      2.86283, 1.85412, 1.04581)))
  expect_that(round(Bout2$sigma2, 5),
    is_equivalent_to(c(1.07199, 0.70925, 0.80124, 0.50868, 0.67254, 0.50322,
      2.75713, 3.05796, 2.42955)))
  expect_that(round(attr(Bout2, "Rhat"), 5),
    is_equivalent_to(c(3.65421,  1.31099,  2.16273,  1.78500,  8.11593)))
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
    is_equivalent_to(c(0.27273, 0.09091, 0.09091, 0.18182, 0.09091,
      0.09091, 0.27273, 0.18182, 0.18182, 0.09091, 0.09091, 0.09091)))
  expect_that(round(pow2[, 2], 5),
    is_equivalent_to(c(0.04645, 0.00000, 0.00000, 0.00729, 0.00000,
      0.00000, 0.04645, 0.00729, 0.00729, 0.00000, 0.00000, 0.00000)))
  expect_that(round(pow2[, 3], 5),
    is_equivalent_to(c(0.52243, 0.25887, 0.25887, 0.39781, 0.25887,
      0.25887, 0.52243, 0.39781, 0.39781, 0.25887, 0.25887, 0.25887)))
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
    is_equivalent_to(c(4.66856, 4.42793, 4.28281, 4.81650, 5.74952,
      4.25399, 4.58407, 5.26500, 4.61872)))
  expect_that(round(Bout2a$mu2, 5),
    is_equivalent_to(c(2.65538, 2.64737, 3.03047, 3.51627, 3.41712,
      3.15124, 3.41061, 3.41331, 3.20443)))
  expect_that(round(Bout2a$nu, 5),
    is_equivalent_to(c(6.93609,  8.59607,  1.14347, 10.86202, 11.87991,
      8.33457,  5.79358,  6.03924,  4.72087)))
  expect_that(round(Bout2a$sigma1, 5),
    is_equivalent_to(c(0.40205, 0.42936, 0.39790, 0.59238, 1.05655,
      2.24298, 2.02154, 1.73155, 1.81695)))
  expect_that(round(Bout2a$sigma2, 5),
    is_equivalent_to(c(0.71316, 0.77907, 0.77214, 0.62478, 0.51189,
      0.83907, 0.65247, 0.53180, 0.94352)))
  expect_that(round(attr(Bout2a, "Rhat"), 5),
    is_equivalent_to(c(1.11267,  2.83848, 1.79341,  2.26465,  0.96004)))
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
    is_equivalent_to(c(1.19359, 1.18673, 1.12595, 1.15407, 0.96022,
      1.39089, 1.31197, 1.63615, 1.65069)))
  expect_that(round(Bout1$nu, 5),
    is_equivalent_to(c(22.76624, 50.06082, 40.31628, 28.49252, 19.36419,
      15.38055,  1.88658,  2.93115, 6.11908)))
  expect_that(round(Bout1$sigma, 5),
    is_equivalent_to(c(0.30475, 0.46296, 0.57263, 0.47637, 0.55622,
      0.77982, 0.48053, 0.30616, 0.33081)))
  expect_that(round(attr(Bout1, "Rhat"), 5),
    is_equivalent_to(c(2.01126, 2.93018, 1.52859)))
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
    is_equivalent_to(c(0.63636, 0.09091, 0.09091, 0.54545, 0.09091,
      0.09091, 0.18182, 0.72727, 0.18182, 0.09091, 0.09091, 0.18182)))
  expect_that(round(pow1[, 2], 5),
    is_equivalent_to(c(0.36680, 0.00000, 0.00000, 0.26839, 0.00000,
      0.00000, 0.00729, 0.47757, 0.00729, 0.00000, 0.00000, 0.00729)))
  expect_that(round(pow1[, 3], 5),
    is_equivalent_to(c(0.89322, 0.25887, 0.25887, 0.81845, 0.25887,
      0.25887, 0.39781, 0.95355, 0.39781, 0.25887, 0.25887, 0.39781)))
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
    is_equivalent_to(c(1.63622, 1.64527, 1.30575, 1.06781, 1.09037,
      1.17121, 1.62116, 1.64793, 1.59050)))
  expect_that(round(Bout1a$nu, 5),
    is_equivalent_to(c(17.35507, 23.91049, 20.09931,  5.07647,  7.55075,
      1.35231,  4.98792,  4.35384, 4.74761)))
  expect_that(round(Bout1a$sigma, 5),
    is_equivalent_to(c(0.35104, 0.34143, 0.46024, 0.56157, 0.38128,
      0.44784, 0.30305, 0.54008, 0.41961)))
  expect_that(round(attr(Bout1a, "Rhat"), 5),
    is_equivalent_to(c(3.51877, 5.22756, 0.98897)))
  expect_that(attr(Bout1a, "n.eff"),
    is_equivalent_to(c(9,  9,  9)))
  expect_that(names(attr(Bout1a, "priors")), equals(c("muM", "muSD", "sigmaMode",
    "sigmaSD", "nuMean", "nuSD")))
})
