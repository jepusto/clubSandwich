library(dplyr)
library(plm)
library(clubSandwich)
rm(list=ls())

# load and clean data

data("AchievementAwardsRCT")

AA_RCT_females <- 
  AchievementAwardsRCT %>%
  filter(sex=="Girl" & year != "1999") %>%
  select(-sex) %>%
  mutate(sibs_4 = siblings >= 4,
         treated2001 = treated * (year=="2001"))

#-------------------------
# ATE model
#-------------------------

ATE_mod <- plm(Bagrut_status ~ year * school_type + 
               father_ed + mother_ed + immigrant + sibs_4 + 
               qrtl + treated2001:half, 
               data = AA_RCT_females, 
               index = c("school_id","student_id"), effect = "individual")

trt_effects <- grepl("treated2001", names(coef(ATE_mod)))

ATE_constraints <- list("ATE - upper half (q = 1)" = which(trt_effects)[2], 
                    "ATE - joint (q = 2)" = which(trt_effects))

# standard CRVE (CR1)
ATE_CR1 <- Wald_test(ATE_mod, constraints = ATE_constraints, 
                     vcov = "CR1", test = "Naive-F")

# AHT tests (CR2)
ATE_CR2 <- Wald_test(ATE_mod, constraints = ATE_constraints, 
                     vcov = "CR2", test = "HTZ")

#---------------------------------------------
# Model with school-by-sector interactions
#---------------------------------------------

school_mod <- plm(Bagrut_status ~ year * school_type + 
                  father_ed + mother_ed + immigrant + sibs_4 + 
                  qrtl + treated2001:half + treated2001:half:school_type, 
                  data = AA_RCT_females, 
                  index = c("school_id","student_id"), effect = "individual")

school_effects <- grepl("school_type.*treated2001", names(coef(school_mod)))

school_constraints <- list("Moderation - upper half (q = 2)" = which(school_effects)[3:4], 
                        "Moderation - joint (q = 3)" = which(school_effects))

# standard CRVE (CR1)
mod_CR1 <- Wald_test(school_mod, constraints = school_constraints, 
                     vcov = "CR1", test = "Naive-F")

# AHT tests (CR2)
mod_CR2 <- Wald_test(school_mod, constraints = school_constraints, 
                     vcov = "CR2", test = "HTZ")

#---------------------------------------
# Arrange test results in a table
#---------------------------------------
ATE_CR1 <- bind_rows(ATE_CR1, .id = "Hypothesis") %>% as.data.frame()
ATE_CR2 <- bind_rows(ATE_CR2, .id = "Hypothesis") %>% as.data.frame()
ATEs <- bind_rows("Standard" = ATE_CR1, "AHT" = ATE_CR2, .id = "Test") %>% 
  arrange(desc(Hypothesis))

mod_CR1 <- bind_rows(mod_CR1, .id = "Hypothesis") %>% as.data.frame()
mod_CR2 <- bind_rows(mod_CR2, .id = "Hypothesis") %>% as.data.frame()
mods <- bind_rows("Standard" = mod_CR1, "AHT" = mod_CR2, .id = "Test") %>% 
  arrange(desc(Hypothesis))

AL_results <- 
  bind_rows(ATEs, mods) %>%
  select(Hypothesis, Test, F = Fstat, df, p = p_val) %>%
  mutate(Hypothesis = ifelse(Test=="AHT",NA,Hypothesis))
