library(AER)
data("CigarettesSW")
Cigs <- within(CigarettesSW, {
  packs_log <- log(packs)
  rprice_log <- log(price/cpi)
  rincome_log <- log(income/population/cpi)
  tdiff <- (taxs - tax)/cpi
})

iv_fit <- ivreg(packs_log ~ rprice_log + rincome_log | 
                  rincome_log + tdiff, data = Cigs)
vcovCR(iv_fit, cluster = Cigs$state, type = "CR2")
coef_test(iv_fit, vcov = "CR2", cluster = Cigs$state)

sur_fit <- lm(cbind(packs_log, rprice_log) ~ rincome_log + tdiff, data = Cigs)
(V_sur <- vcovCR(sur_fit, cluster = Cigs$state, type = "CR2"))
coef_test(sur_fit, vcov = "CR2", cluster = Cigs$state)

(a <- coef_CS(sur_fit)[["rprice_log:tdiff"]])
(b <- coef_CS(sur_fit)[["packs_log:tdiff"]])
(d <- coef_CS(iv_fit)[["rprice_log"]])
b / a
all.equal(d, b / a)

(V_iv <- vcovCR(iv_fit, cluster = Cigs$state, type = "CR2")["rprice_log","rprice_log"])
(V_delta <- (V_sur["packs_log:tdiff","packs_log:tdiff"] + 
              d^2 * V_sur["rprice_log:tdiff","rprice_log:tdiff"] - 
              2 * d * V_sur["rprice_log:tdiff","packs_log:tdiff"]) / a^2)
all.equal(V_iv, V_delta)
