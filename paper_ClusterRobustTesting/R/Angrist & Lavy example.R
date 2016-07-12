library(dplyr)
library(tidyr)
library(foreign)
library(plm)
library(clubSandwich)

#-----------------------------------------
# read in and munge data
#-----------------------------------------
# setwd("paper_ClusterRobustTesting")

filenames <- paste0("data/AngristLavy_AERdata/base",c("99","00","01","02"), ".dta")
panel_dat_list <- lapply(filenames, read.dta)
panel_dat <- bind_rows(panel_dat_list)

panel_dat <- within(panel_dat, {
  yr <- factor(year, levels = c(99,0,1,2), labels = 1999:2002)
  school_type <- factor(ifelse(semrel, "Religious", ifelse(semarab, "Arab","Secular")))
  student_id <- paste0(yr, "-", student_id)
  treated2001 <- (yr=="2001") * treated
  sex <- ifelse(boy==1, "Boy","Girl")
  ah4 <- m_ahim >= 4
})

panel_dat <- cbind(panel_dat, model.matrix(~ 0 + yr, data = panel_dat))

qntl <- function(x, groups = 4) {
  qtl <- quantile(x, (0:groups) / groups)
  cut(x, qtl, labels = 1:groups, include.lowest = TRUE, right = FALSE, ordered_result = TRUE)
}

filter(panel_dat, yr!=1999) %>%
  group_by(yr, sex) %>%
  mutate(qrtl = qntl(lagscore, 4),
         half = qntl(lagscore, 2)) %>%
  ungroup() %>%
  filter(sex=="Girl") %>%
  select(-yr1999, -boy, -sex) %>%
  arrange(school_id, student_id) ->
  panel_dat

panel_fit <- plm(zakaibag ~ yr2000:school_type + yr2001:school_type + 
                   educav + educem + ole5 + ah4 + qrtl + treated2001:half, 
                 data = panel_dat, index = c("school_id","student_id"), effect = "individual")
panel_constraints <<- grepl("treated2001", names(coef(panel_fit)))

# fit model after absorbing year effects
X1_mat <- model.matrix(panel_fit)
R1_coefs <- c(1:7,14:15)
R1_mat <- residuals(lm(X1_mat[,R1_coefs] ~ X1_mat[,-R1_coefs]))
absorb_fit1 <- lm(panel_dat$zakaibag ~ 0 + R1_mat)
# cbind(coef(absorb_fit1), coef(panel_fit)[R1_coefs])
# all.equal(coef(absorb_fit1), coef(panel_fit)[R1_coefs], check.attributes = FALSE)

# test ATE by half
t_CR1 <- coef_test(panel_fit, vcov="CR1", test = "naive-t")[panel_constraints,]
t_CR2 <- coef_test(panel_fit, vcov="CR2", test = "Satterthwaite")[panel_constraints,]
t_CR2A <- coef_test(absorb_fit1, vcov="CR2", cluster = panel_dat$school_id, test = "Satterthwaite")[8:9,]

joint_CR1 <- Wald_test(panel_fit, constraints = panel_constraints, vcov = "CR1", test = "Naive-F")
joint_CR2 <- Wald_test(panel_fit, constraints = panel_constraints, vcov = "CR2", test = "HTZ")
joint_CR2A <- Wald_test(absorb_fit1, constraints = 8:9, vcov = "CR2", cluster = panel_dat$school_id, test = "HTZ")

panel_school <- plm(zakaibag ~ yr2000:school_type + yr2001:school_type + educav + educem + ole5 + ah4 + qrtl 
                    + treated2001:half + treated2001:half:school_type, 
                    data = panel_dat, index = c("school_id","student_id"), effect = "individual")
school_constraints <- grepl("school_type.*treated2001", names(coef(panel_school)))

# fit model after absorbing year effects
X2_mat <- model.matrix(panel_school)
R2_coefs <- c(1:7,14:19)
R2_mat <- residuals(lm(X2_mat[,R2_coefs] ~ X2_mat[,-R2_coefs]))
absorb_fit2 <- lm(panel_dat$zakaibag ~ 0 + R2_mat)
# cbind(coef(absorb_fit2), coef(panel_school)[R2_coefs])
# all.equal(coef(absorb_fit2), coef(panel_school)[R2_coefs], check.attributes = FALSE)

# test moderation
constraint_list <- list(lower = 16:17, upper = 18:19, joint = 16:19)
mod_CR1 <- Wald_test(panel_school, constraint_list, vcov = "CR1", test = "Naive-F")
mod_CR2 <- Wald_test(panel_school, constraint_list, vcov = "CR2", test = "HTZ")
mod_CR2A <- Wald_test(absorb_fit2, constraints = list(lower = 10:11, upper = 12:13, joint = 10:13), 
                     vcov = "CR2", cluster = panel_dat$school_id, test = "HTZ")

t_CR1 <- within(t_CR1, {
    df <- 34
    hypothesis <- rownames(t_CR1)
    CR <- "CR1"
    Test <- "ad hoc"
    p_val <- p_t
  })
t_CR2 <- within(t_CR2, {
  hypothesis <- rownames(t_CR2)
  CR <- "CR2"
  Test <- "AHT"
  p_val <- p_Satt
})
t_CR2A <- within(t_CR2A, {
  hypothesis <- rownames(t_CR2)
  CR <- "CR2"
  Test <- "AHT*"
  p_val <- p_Satt
})

t_tests <- 
  bind_rows(select(t_CR1, hypothesis, CR, Test, beta, SE, df, p_val),
            select(t_CR2, hypothesis, CR, Test, beta, SE, df, p_val), 
            select(t_CR2A, hypothesis, CR, Test, beta, SE, df, p_val)) %>%
  mutate(Fstat = (beta / SE)^2) %>%
  filter(hypothesis=="treated2001:half2") %>%
  mutate(Hypothesis = "ATE - upper half (q = 1)") %>%
  select(Hypothesis, Test, Fstat, df, p = p_val)

joint_tests <- 
  bind_rows("ad hoc" = joint_CR1, 
            "AHT" = joint_CR2, 
            "AHT*" = joint_CR2A, .id = "Test") %>%
  mutate(Hypothesis = "ATE - joint (q = 2)") %>%
  select(Hypothesis, Test, Fstat, df, p = p_val)


mod_tests <- 
  bind_rows("ad hoc" = bind_rows(mod_CR1, .id = "hypothesis"),
            "AHT" = bind_rows(mod_CR2, .id = "hypothesis"),
            "AHT*" = bind_rows(mod_CR2A, .id = "hypothesis"), .id = "Test") %>%
  filter(hypothesis %in% c("upper","joint")) %>%
  arrange(desc(hypothesis), Test) %>%
  mutate(Hypothesis = ifelse(hypothesis=="upper", "Moderation - upper half (q = 2)", "Moderation - joint (q = 4)")) %>%
  select(Hypothesis, Test, Fstat, df, p = p_val)

AL_results <- 
  bind_rows(t_tests, joint_tests, mod_tests) %>%
  mutate(df = round(df, 2),
         p = round(p, 5))

AL_results$Hypothesis[-seq(1,12,3)] <- NA
