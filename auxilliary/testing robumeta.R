library(robumeta)
library(clubSandwich)
data(dropoutPrevention)

dp_subset <- subset(dropoutPrevention, big_study==TRUE)

m3_corr_small <- robu(LOR1 ~ study_design + attrition + group_equivalence + adjusted
                      + outcome + evaluator_independence
                      + male_pct + white_pct + average_age
                      + implementation_quality + program_site + duration + service_hrs, 
                      data = dp_subset, studynum = studyID, var.eff.size = varLOR, modelweights = "CORR")

m3_corr_small
coef_test(m3_corr_small, vcov = "CR2")


m3_hier_small <- robu(LOR1 ~ study_design + attrition + group_equivalence + adjusted
                + outcome + evaluator_independence
                + male_pct + white_pct + average_age
                + implementation_quality + program_site + duration + service_hrs, 
                data = dp_subset, studynum = studyID, var.eff.size = varLOR, modelweights = "HIER")

m3_hier_small

ord <- order(order(m3_hier_small$study_orig_id))

dp_subset$user_wts <- m3_hier_small$data.full$r.weights[ord]

m3_user_small <- robu(LOR1 ~ study_design + attrition + group_equivalence + adjusted
                      + outcome + evaluator_independence
                      + male_pct + white_pct + average_age
                      + implementation_quality + program_site + duration + service_hrs, 
                      data = dp_subset, studynum = studyID, var.eff.size = varLOR, userweights = user_wts)

m3_user_small
coef_test(m3_user_small, vcov = "CR2")
