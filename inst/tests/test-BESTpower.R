
# Tests for BESTpower



context("BESTpower_Pro")

test_that("BESTpower with 2 groups gives same output",  {
  proData <- makeData(mu1=108, sd1=17, mu2=100, sd2=15, nPerGrp=20, 
                         pcntOut=10, sdOutMult=2.0, rnd.seed=1,
                         showPlot=FALSE)
  proMCMC <- BESTmcmc(proData$y1, proData$y2, numSavedSteps=9,
      burnInSteps = 1, verbose=FALSE, rnd.seed=2)  
  pow2 <- BESTpower(proMCMC, N1=10, N2=10,
               ROPEm=c(-2,2) , ROPEsd=c(-2,2) , ROPEeff=c(-0.5,0.5) , 
               maxHDIWm=25.0 , maxHDIWsd=10.0 , maxHDIWeff=1.0 ,
               nRep=9, mcmcLength=1000, verbose=0, rnd.seed=3) 
                       
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
    is_equivalent_to(c(0.18182, 0.09091, 0.09091, 0.09091, 0.09091,
      0.09091, 0.09091, 0.09091, 0.09091, 0.09091, 0.09091, 0.09091)))
  expect_that(round(pow2[, 2], 5), 
    is_equivalent_to(c(0.00729, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)))
  expect_that(round(pow2[, 3], 5), 
    is_equivalent_to(c(0.39781, 0.25887, 0.25887, 0.25887, 0.25887,
      0.25887, 0.25887, 0.25887, 0.25887, 0.25887, 0.25887, 0.25887)))
})

test_that("BESTpower with 1 group gives same output",  {
  proData <- makeData(mu1=108, sd1=17, nPerGrp=20, 
                         pcntOut=10, sdOutMult=2.0, rnd.seed=4,
                         showPlot=FALSE)
  proMCMC <- BESTmcmc(proData$y1, proData$y2, numSavedSteps=9,
      burnInSteps = 1, verbose=FALSE, rnd.seed=2)  
  pow1 <- BESTpower(proMCMC, N1=10, N2=10,
               ROPEm=c(-2,2) , ROPEsd=c(-2,2) , ROPEeff=c(-0.5,0.5) , 
               maxHDIWm=25.0 , maxHDIWsd=10.0 , maxHDIWeff=1.0 ,
               nRep=9, mcmcLength=1000, verbose=0, rnd.seed=3) 
                       
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
    is_equivalent_to(c(0.90909, 0.09091, 0.09091, 0.36364, 0.90909,
      0.09091, 0.09091, 0.09091, 0.90909, 0.09091, 0.09091, 0.09091)))
  expect_that(round(pow1[, 2], 5), 
    is_equivalent_to(c(0.74113, 0, 0, 0.10678, 0.74113, 0, 0, 0, 0.74113, 0, 0, 0)))
  expect_that(round(pow1[, 3], 5), 
    is_equivalent_to(c(1, 0.25887, 0.25887, 0.6332, 1, 0.25887, 0.25887,
      0.25887, 1, 0.25887, 0.25887, 0.25887)))
})

