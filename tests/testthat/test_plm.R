# context("plm objects")
# 
# library(plm)
# 
# data("Produc", package = "plm")
# Produc$cluster <- sample(LETTERS[1:10], size = nrow(Produc), replace=TRUE)
# Produc_scramble <- Produc[sample(nrow(Produc)),]
# 
# test_that("individual effects agree with lm", {
#   plm_individual <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
#                              data = Produc_scramble, index = c("state","year"), 
#                              effect = "individual", model = "within")
#   lm_individual <- lm(log(gsp) ~ 0 + state + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
#   coef_names <- names(coef(plm_individual)) 
#   coef_index <- names(coef(lm_individual)) %in% coef_names
#   
#   expect_equal(vcovCR(plm_individual, type="CR0")[coef_names,coef_names], 
#                vcovCR(lm_individual, cluster = Produc$state, type = "CR0")[coef_index,coef_index])
#   expect_equal(vcovCR(plm_individual, type="CR1")[coef_names,coef_names], 
#                vcovCR(lm_individual, cluster = Produc$state, type = "CR1")[coef_index,coef_index])
#   expect_equal(vcovCR(plm_individual, type="CR2")[coef_names,coef_names], 
#                vcovCR(lm_individual, cluster = Produc$state, type = "CR2")[coef_index,coef_index])
# })
# 
# test_that("time effects agree with lm", {
#   plm_time <- plm::plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
#                        data = Produc_scramble, index = c("state","year"), 
#                        effect = "time", model = "within")
#   lm_time <- lm(log(gsp) ~ 0 + factor(year) + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
#   coef_names <- names(coef(plm_time)) 
#   coef_index <- names(coef(lm_time)) %in% coef_names
#   
#   expect_equal(vcovCR(plm_time, type="CR0")[coef_names,coef_names], 
#                vcovCR(lm_time, cluster = Produc$year, type = "CR0")[coef_index,coef_index])
#   expect_equal(vcovCR(plm_time, type="CR1")[coef_names,coef_names], 
#                vcovCR(lm_time, cluster = Produc$year, type = "CR1")[coef_index,coef_index])
#   expect_equal(vcovCR(plm_time, type="CR2")[coef_names,coef_names], 
#                vcovCR(lm_time, cluster = Produc$year, type = "CR2")[coef_index,coef_index])
# })
# 
# test_that("two-way effects agree with lm", {
#   plm_twoways <- plm::plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
#                           data = Produc_scramble, index = c("state","year"), 
#                           effect = "twoways", model = "within")
#   lm_twoways <- lm(log(gsp) ~ 0 + state + factor(year) + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
#   coef_names <- names(coef(plm_twoways)) 
#   coef_index <- names(coef(lm_twoways)) %in% coef_names
#   
#   # clustering on individual
#   expect_equal(vcovCR(plm_twoways, cluster = "individual", type="CR0")[coef_names,coef_names], 
#                vcovCR(lm_twoways, cluster = Produc$state, type = "CR0")[coef_index,coef_index])
#   expect_equal(vcovCR(plm_twoways, cluster = "individual", type="CR1")[coef_names,coef_names], 
#                vcovCR(lm_twoways, cluster = Produc$state, type = "CR1")[coef_index,coef_index])
#   expect_equal(vcovCR(plm_twoways, cluster = "individual", type="CR2")[coef_names,coef_names], 
#                vcovCR(lm_twoways, cluster = Produc$state, type = "CR2")[coef_index,coef_index])
# 
#   # clustering on time
#   expect_equal(vcovCR(plm_twoways, cluster = "time", type="CR0")[coef_names,coef_names], 
#                vcovCR(lm_twoways, cluster = Produc$year, type = "CR0")[coef_index,coef_index])
#   expect_equal(vcovCR(plm_twoways, cluster = "time", type="CR1")[coef_names,coef_names], 
#                vcovCR(lm_twoways, cluster = Produc$year, type = "CR1")[coef_index,coef_index])
#   expect_equal(vcovCR(plm_twoways, cluster = "time", type="CR2")[coef_names,coef_names], 
#                vcovCR(lm_twoways, cluster = Produc$year, type = "CR2")[coef_index,coef_index])
# 
#   # clustering on a randomly generated factor
#   expect_equal(vcovCR(plm_twoways, cluster = Produc_scramble$cluster, type="CR0")[coef_names,coef_names], 
#                vcovCR(lm_twoways, cluster = Produc$cluster, type = "CR0")[coef_index,coef_index])
#   expect_equal(vcovCR(plm_twoways, cluster = Produc_scramble$cluster, type="CR1")[coef_names,coef_names], 
#                vcovCR(lm_twoways, cluster = Produc$cluster, type = "CR1")[coef_index,coef_index])
#   expect_equal(vcovCR(plm_twoways, cluster = Produc_scramble$cluster, type="CR2")[coef_names,coef_names], 
#                vcovCR(lm_twoways, cluster = Produc$cluster, type = "CR2")[coef_index,coef_index])
#   
# })

# test for equality with HC when cluster = rownames(data)
# test cluster specification
# test target matrix specification
# test inverse_var detection