
data(STAR, package = "AER")

library(clubSandwich)
library(dplyr)
library(tidyr)

# pull out first grade information

STAR_1st_urban <- 
  STAR %>%
  select(gender, ethnicity, birth, ends_with("1")) %>%
  filter(!is.na(star1)) %>%
  mutate(student_id = row_number(), 
         schoolid1 = factor(schoolid1)) %>%
  
# limit to urban schools
  filter(school1=="urban") %>%
  
# re-shape to long format
  gather(key = "test", value = "score", read1, math1)

# seemingly unrelated regression 

lm_STAR <- lm(score ~ test:(0 + schoolid1 + gender + ethnicity + birth + lunch1 + star1), data = STAR_1st_urban)

# joint test for small class effects on reading and math

small_class_coefs <- which(grepl("star1small", names(coef(lm_STAR))))


# clustering by student 

Wald_test(lm_STAR, constraints = small_class_coefs, 
          vcov = "CR2", cluster = STAR_1st_urban$student_id, 
          test = "HTZ")

# clustering by school 

Wald_test(lm_STAR, constraints = small_class_coefs, 
          vcov = "CR2", cluster = STAR_1st_urban$schoolid1, 
          test = "HTZ")

