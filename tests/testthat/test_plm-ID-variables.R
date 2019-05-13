context("plm objects - ID variables")
set.seed(20190513)

library(plm, quietly=TRUE)

data("Produc", package = "plm")
Produc <- Produc[sample(nrow(Produc)),]
Produc$cluster <- sample(LETTERS[1:10], size = nrow(Produc), replace=TRUE)


plm_individual <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                      data = Produc, index = "state",
                      effect = "individual", model = "within")

lm_individual <- lm(log(gsp) ~ 0 + state + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
individual_names <- names(coef(plm_individual)) 
individual_index <- names(coef(lm_individual)) %in% individual_names

lm_CR0 <- vcovCR(lm_individual, cluster = Produc$state, type = "CR0")[individual_index,individual_index]
lm_CR1 <- vcovCR(lm_individual, cluster = Produc$state, type = "CR1")[individual_index,individual_index]
lm_CR2 <- vcovCR(lm_individual, cluster = Produc$state, type = "CR2")[individual_index,individual_index]

plm_CR0 <- vcovCR(plm_individual, type="CR0")[individual_names,individual_names]

test_that("individual effects agree with lm under automatic clustering", {
  plm_CR1 <- vcovCR(plm_individual, type="CR1")[individual_names,individual_names]
  plm_CR2 <- vcovCR(plm_individual, type="CR2")[individual_names,individual_names]
  
  expect_equal(plm_CR0, lm_CR0)
  expect_equal(plm_CR1, lm_CR1)
  expect_equal(plm_CR2, lm_CR2)
})

test_that("individual effects agree with lm under explicit clustering", {
  plm_CR1 <- vcovCR(plm_individual, cluster = Produc$state, type="CR1")[individual_names,individual_names]
  plm_CR2 <- vcovCR(plm_individual, cluster = Produc$state, type="CR2")[individual_names,individual_names]

  expect_equal(plm_CR1, lm_CR1)
  expect_equal(plm_CR2, lm_CR2)
})

test_that("individual effects agree with lm under random clustering", {
  
  lm_CR1 <- vcovCR(lm_individual, cluster = Produc$cluster, type = "CR1")[individual_index,individual_index]
  lm_CR2 <- vcovCR(lm_individual, cluster = Produc$cluster, type = "CR2")[individual_index,individual_index]

  plm_CR1 <- vcovCR(plm_individual, cluster = Produc$cluster, type="CR1")[individual_names,individual_names]
  plm_CR2 <- vcovCR(plm_individual, cluster = Produc$cluster, type="CR2")[individual_names,individual_names]
  
  expect_equal(plm_CR1, lm_CR1)
  expect_equal(plm_CR2, lm_CR2)
})

test_that("CR0 and CR1S agree with arellano vcov", {
  expect_equal(vcovHC(plm_individual, method="arellano", type = "HC0", cluster = "group"), 
               as.matrix(plm_CR0), check.attributes = FALSE)
  expect_equal(vcovHC(plm_individual, method="arellano", type = "sss", cluster = "group"), 
               as.matrix(vcovCR(plm_individual, type = "CR1S")), check.attributes = FALSE)
})



test_that("plm works for Yuki Takahashi's reprex.",{
  N <- 100
  id <- rep(1:N, 2)
  gid <- rep(1:(N/2), 4)
  Trt <- rep(c(0,1), each = N)
  a <- rep(rnorm(N, mean=0, sd=0.005), 2)
  gp <- rep(rnorm(N / 2, mean=0, sd=0.0005), 4)
  u <- rnorm(N * 2, mean=0, sd=0.05)
  Ylatent <- -0.05 * Trt + gp + a + u
  
  Data <- data.frame(
    Y = ifelse(Ylatent > 0, 1, 0),
    id, gid, Trt
  )
  
  fe_fit <- plm(formula = Y ~ Trt, data = Data, 
                model = "within", index = "id", 
                effect = "individual", 
                singular.ok = FALSE)
  
  
  implicit <- vcovCR(fe_fit, type = "CR2")
  explicit <- vcovCR(fe_fit, cluster=Data$id, type = "CR2") 

  expect_equal(implicit, explicit)
  expect_s3_class(vcovCR(fe_fit, cluster=Data$gid, type = "CR2"), "vcovCR") 
  
})



test_that("Clustering works for various ways of specifying unit and time indices in plm.", {
  
  data("Grunfeld", package = "plm")
  Grunfeld$cluster <- sample(LETTERS[1:10], size = nrow(Grunfeld), replace=TRUE)
  rearrange <- mapply(function(s,b) seq(s, nrow(Grunfeld), b), 1:10, 11:20)
  rearrange <- unique(unlist(rearrange))
  rearrange <- c(rearrange, setdiff(1:200, rearrange))
  Grunfeld_scramble <- Grunfeld[rearrange,]
  Grunfeld_pdata <- pdata.frame(Grunfeld_scramble, index = c("firm","year"))
  
  plm_pdata <- plm(inv ~ value + capital, data = Grunfeld_pdata, model="within")
  plm_numeric <- plm(inv ~ value + capital, data = Grunfeld, index = 10, model="within")
  plm_noindex <- plm(inv ~ value + capital, data = Grunfeld_scramble, model="within")
  plm_oneindex <- plm(inv ~ value + capital, data = Grunfeld_scramble, index = "firm", model="within")
  plm_twoindex <- plm(inv ~ value + capital, data = Grunfeld_scramble, index = c("firm","year"), model="within")

  CR_types <- paste0("CR",0:3)
  
  # auto clustering
  vcovCRs <- function(model, types) lapply(types, function(x) vcovCR(model, type = x))
  CR_pdata <- vcovCRs(plm_pdata, CR_types)
  expect_equivalent(CR_pdata, vcovCRs(plm_numeric, CR_types))
  expect_equivalent(CR_pdata, vcovCRs(plm_noindex, CR_types))
  expect_equivalent(CR_pdata, vcovCRs(plm_oneindex, CR_types))
  expect_equivalent(CR_pdata, vcovCRs(plm_twoindex, CR_types))
  
  # manual clustering on firm
  vcovCRs <- function(model, types, cluster) lapply(types, function(x) vcovCR(model, type = x, cluster = cluster))
  expect_equivalent(CR_pdata, vcovCRs(plm_numeric, CR_types, cluster = Grunfeld$firm))
  expect_equivalent(CR_pdata, vcovCRs(plm_noindex, CR_types, cluster = Grunfeld_scramble$firm))
  expect_equivalent(CR_pdata, vcovCRs(plm_oneindex, CR_types, cluster = Grunfeld_scramble$firm))
  expect_equivalent(CR_pdata, vcovCRs(plm_twoindex, CR_types, cluster = Grunfeld_scramble$firm))

  # manual clustering on arbitrary id
  CR_pdata <- vcovCRs(plm_pdata, CR_types, cluster = Grunfeld_pdata$cluster)
  expect_equivalent(CR_pdata, vcovCRs(plm_numeric, CR_types, cluster = Grunfeld$cluster))
  expect_equivalent(CR_pdata, vcovCRs(plm_noindex, CR_types, cluster = Grunfeld_scramble$cluster))
  expect_equivalent(CR_pdata, vcovCRs(plm_oneindex, CR_types, cluster = Grunfeld_scramble$cluster))
  expect_equivalent(CR_pdata, vcovCRs(plm_twoindex, CR_types, cluster = Grunfeld_scramble$cluster))
  
})

