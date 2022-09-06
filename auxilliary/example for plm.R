library(plm)
data("Produc")

plm_individual <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                      data = Produc, index = c("state","year","region"), 
                      effect = "individual", model = "random")

plm_vcov <- vcovHC(plm_individual, method="arellano", type = "HC0", cluster = "group")
club_vcov <- vcovCR(plm_individual, cluster = "individual", type = "CR0")
all.equal(plm_vcov, as.matrix(club_vcov), check.attributes = FALSE)

plm_nested <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                  data = Produc, index = c("state","year","region"), 
                  effect = "nested", model = "random")

plm_vcov <- vcovHC(plm_nested, method="arellano", type = "HC0", cluster = "group")
club_vcov <- vcovCR(plm_nested, cluster = "group", type = "CR0")
all.equal(plm_vcov, as.matrix(club_vcov), check.attributes = FALSE)
