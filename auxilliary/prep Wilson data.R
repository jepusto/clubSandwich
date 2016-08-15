library(plyr)
library(foreign)
library(robumeta)
rm(list=ls())

Wilson <- droplevels(read.spss("auxilliary/Dropout Dataset for Tipton and Pustejovsky.sav", to.data.frame = TRUE))
Wilson <- rename(Wilson, replace = 
                   c("attritionimp" = "attrition", "LGEimp" = "group_equivalence", 
                     "es50r" = "adjusted","g10i" = "evaluator_independence",
                     "g42_TXimp" = "male_pct", "g43a_TXimp" = "white_pct", 
                     "g46_TXimp" = "average_age", "g30" = "implementation_quality", 
                     "g20imp" = "duration", "g24imp" = "service_hrs",
                     "g6r" = "Program_type"))

Wilson <- within(Wilson, {
  levels(implementation_quality) <- c("Clear problems","Possible problems","No apparent problems")
  levels(evaluator_independence) <- c("Independent","Indirect, influential","Planning","Delivery")
  program_site <- factor(ifelse(class == 1, "school classroom", 
                               ifelse(outclass == 1, "school, outside of classroom",
                                      ifelse(mixed == 1, "mixed", "community"))))
  outcome <- factor(ifelse(dropout == 1 , "dropout", 
                        ifelse(grad == 1, "graduation", 
                               ifelse(gradged == 1, "graduation/ged", "enrolled"))))
  study_design <- factor(ifelse(random==1, "Randomized", 
                             ifelse(matched == 1, "Matched", "Non-random, non-matched")))
  studyID <- factor(studyid)
  studySample <- factor(StudyGroup)
})

Wilson <- arrange(Wilson, studyID, studySample)
Wilson <- subset(Wilson, select = c(LOR1, varLOR, studyID, studySample,
         study_design, outcome, evaluator_independence, implementation_quality, 
         program_site, attrition, group_equivalence, adjusted, 
         male_pct, white_pct, average_age, duration, service_hrs))

ES_per_sample <- ddply(Wilson, .(studyID, studySample), nrow)
table(table(ES_per_sample$studyID))
ES_per_study <- ddply(Wilson, .(studyID), nrow)
table(ES_per_study$V1)
big_studies <- subset(ES_per_study, V1 > 2)$studyID
Wilson$big_study <- Wilson$studyID %in% big_studies

dropoutPrevention <- Wilson
save(dropoutPrevention, file = "data/dropoutPrevention.RData", compress = "xz")
write.csv(dropoutPrevention, file = "auxilliary/dropoutPrevention.csv")

names(dropoutPrevention)
table(dropoutPrevention$study_design)
table(dropoutPrevention$outcome)
table(dropoutPrevention$evaluator_independence)
table(dropoutPrevention$implementation_quality)
table(dropoutPrevention$program_format)
summary(dropoutPrevention$attrition)
summary(dropoutPrevention$group_equivalence)
table(dropoutPrevention$adjusted)
summary(dropoutPrevention$male_pct)
summary(dropoutPrevention$white_pct)
summary(dropoutPrevention$average_age)
summary(dropoutPrevention$duration)
summary(dropoutPrevention$service_hrs)
table(dropoutPrevention$big_study)


library(robumeta)
library(clubSandwich)
data("dropoutPrevention")

m3_hier_full <- robu(LOR1 ~ study_design + attrition + group_equivalence + adjusted
                      + outcome + evaluator_independence
                      + male_pct + white_pct + average_age
                      + implementation_quality + program_site + duration + service_hrs, 
                      data = dropoutPrevention, studynum = studyID, var.eff.size = varLOR, modelweights = "HIER")
contrast_list <- list("Study design" = 2:3, 
                      "Outcome measure" = 7:9,
                      "Evaluator independence" = 10:12,
                      "Implmentation quality" = 16:17,
                      "Program format" = 18:20)

dropout_full <- Wald_test(m3_hier_full, constraints = contrast_list, 
                           vcov = "CR2", cluster = dropoutPrevention$studyID, 
                          test = "HTZ")


dp_subset <- subset(dropoutPrevention, big_study==TRUE)

m3_hier_subset <- robu(LOR1 ~ study_design + attrition + group_equivalence + adjusted
                        + outcome + evaluator_independence
                        + male_pct + white_pct + average_age
                        + implementation_quality + program_site + duration + service_hrs, 
                        data = dp_subset, studynum = studyID, var.eff.size = varLOR, modelweights = "HIER")

contrast_list <- list("Study design" = 2:3, 
                      "Outcome measure" = 7:9,
                      "Evaluator independence" = 10:11,
                      "Implmentation quality" = 15:16,
                      "Program format" = 17:19)

dropout_subset <- Wald_test(m3_hier_subset, constraints = contrast_list, 
                           vcov = "CR2", cluster = dp_subset$studyID,
                           test = "HTZ")
dropout_subset
