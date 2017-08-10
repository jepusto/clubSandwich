library(clubSandwich)
library(nlme)

set.seed(20170802)

# sample sizes at each level
n_districts <- 10
n_schools_per <- rnbinom(n_districts, size = 4, prob = 0.3)
n_schools <- sum(n_schools_per)
n_students_per <- 10
n_students <- n_schools * n_students_per

# identifiers for each level
district_id <- factor(rep(1:n_districts, n_schools_per * n_students_per))
school_id <- factor(rep(1:sum(n_schools_per), each = n_students_per))
student_id <- 1:n_students

# simulated outcome
Y_ijk <- rnorm(n_districts)[district_id] + rnorm(n_schools)[school_id] + rnorm(n_students)
dat <- data.frame(district_id, school_id, student_id, Y = Y_ijk)
head(dat)

#----------------------------------
# fit three-level model
#----------------------------------

lme_3level <- lme(Y ~ 1, random = ~ 1 | district_id / school_id, data = dat)
summary(lme_3level)

# clustering variable is assumed to be at highest level 
coef_test(lme_3level, vcov = "CR2") 

# specifying the clustering variable manually returns the same result
coef_test(lme_3level, vcov = "CR2", cluster = dat$district_id) 

#----------------------------------
# fit two-level model
#----------------------------------

lme_2level <- lme(Y ~ 1, random = ~ 1 | school_id, data = dat)
summary(lme_2level) # note that SE is too small!

# now clustering variable is assumed to be at school level 
coef_test(lme_2level, vcov = "CR2") # SE is still too small!

# specifying the clustering variable manually returns the same result
coef_test(lme_2level, vcov = "CR2", cluster = dat$school_id) 

# specifying the clustering variable at the district level 
coef_test(lme_2level, vcov = "CR2", cluster = dat$district_id) 
