skip_if_not_installed("carData")

data(Baumann, package = "carData")
Baumann$id <- 1:nrow(Baumann)

test_that("Wald_test words with multiple comparisons adjustment and a single test type (Baumann data).", {
  
  Baumann_fit <- lm(post.test.1 ~ 0 + group + pretest.1, data=Baumann)
  
  Wald5 <- Wald_test(
    Baumann_fit,
    constraints = constrain_pairwise("group", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Baumann$id,
    test = "HTZ"
  )
  
  Wald6 <- Wald_test(
    Baumann_fit,
    constraints = constrain_pairwise("group", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Baumann$id,
    test = "HTZ",
    adjustment_method = "none"
  )
  
  # check that explicitly stating default does not affect functionality
  lapply(Wald5, expect_s3_class, class = "Wald_test_clubSandwich")
  lapply(Wald6, expect_s3_class, class = "Wald_test_clubSandwich")
  expect_equal(Wald5, Wald6) 
  
  # change Wald6 to have hochberg adjustment
  Wald6 <- Wald_test(
    Baumann_fit,
    constraints = constrain_pairwise("group", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Baumann$id,
    test = "HTZ",
    adjustment_method = "hochberg"
  )
  
  Wald5_p_values <- sapply(Wald5, function(x) x$p_val) # extract p-values of Wald5
  Wald6_p_values <- sapply(Wald6, function(x) x$p_val) # extract p-values of Wald6
  
  expect_false(all(Wald5_p_values == Wald6_p_values))
  
  # get adjusted_p_values for Wald5, formatted the same as extracted p-values from Wald6
  Wald5_adjusted_p <- p.adjust(Wald5_p_values, method = "hochberg")
  expect_equal(Wald5_adjusted_p, Wald6_p_values)
  
  # Now using tidy = TRUE
  Wald7 <- Wald_test(
    Baumann_fit,
    constraints = constrain_pairwise("group", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Baumann$id,
    test = "EDF",
    tidy = TRUE
  )
  
  Wald8 <- Wald_test(
    Baumann_fit,
    constraints = constrain_pairwise("group", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Baumann$id,
    test = "EDF",
    adjustment_method = "none",
    tidy = TRUE
  )
  
  expect_equal(Wald7,Wald8)
  
  Wald9 <- Wald_test(
    Baumann_fit,
    constraints = constrain_pairwise("group", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Baumann$id,
    test = "EDF",
    adjustment_method = "holm",
    tidy = TRUE
  )
  
  Wald8_p_adjusted <- tapply(Wald8$p_val, Wald8$test, p.adjust, method = "holm", simplify = FALSE)
  Wald8_adjusted <- Wald8
  Wald8_adjusted$p_val <- unsplit(Wald8_p_adjusted, Wald8$test)
  
  expect_equal(Wald8_adjusted, Wald9)
  
})

test_that("Wald_test words with multiple comparisons adjustment and multiple test types (Baumann data).", {
  
  Baumann_fit <- lm(post.test.1 ~ 0 + group + pretest.1, data=Baumann)
  
  Wald5 <- Wald_test(
    Baumann_fit,
    constraints = constrain_pairwise("group", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Baumann$id,
    test = c("HTZ","chi-sq")
  )
  
  Wald6 <- Wald_test(
    Baumann_fit,
    constraints = constrain_pairwise("group", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Baumann$id,
    test = c("HTZ","chi-sq"),
    adjustment_method = "none"
  )
  
  # check that explicitly stating default does not affect functionality
  lapply(Wald5, expect_s3_class, class = "Wald_test_clubSandwich")
  lapply(Wald6, expect_s3_class, class = "Wald_test_clubSandwich")
  expect_equal(Wald5, Wald6) 
  
  # change Wald6 to have hochberg adjustment
  Wald6 <- Wald_test(
    Baumann_fit,
    constraints = constrain_pairwise("group", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Baumann$id,
    test = c("HTZ","chi-sq"),
    adjustment_method = "hochberg"
  )
  
  Wald5_p_values <- sapply(Wald5, function(x) x$p_val) # extract p-values of Wald5
  Wald6_p_values <- sapply(Wald6, function(x) x$p_val) # extract p-values of Wald6
  
  expect_false(all(Wald5_p_values == Wald6_p_values))
  
  # get adjusted_p_values for Wald5, formatted the same as extracted p-values from Wald6
  Wald5_adjusted_p <- apply(Wald5_p_values, 1, p.adjust, method = "hochberg", simplify = FALSE)
  Wald5_adjusted_p <- do.call(rbind, Wald5_adjusted_p)
  expect_equal(Wald5_adjusted_p, Wald6_p_values)
  
  # Now using tidy = TRUE
  Wald7 <- Wald_test(
    Baumann_fit,
    constraints = constrain_pairwise("group", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Baumann$id,
    test = c("HTA","EDF","EDT"),
    tidy = TRUE
  )
  
  Wald8 <- Wald_test(
    Baumann_fit,
    constraints = constrain_pairwise("group", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Baumann$id,
    test = c("HTA","EDF","EDT"),
    adjustment_method = "none",
    tidy = TRUE
  )
  
  expect_equal(Wald7,Wald8)
  
  Wald9 <- Wald_test(
    Baumann_fit,
    constraints = constrain_pairwise("group", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Baumann$id,
    test = c("HTA","EDF","EDT"),
    adjustment_method = "bonferroni",
    tidy = TRUE
  )
  
  Wald8_p_adjusted <- tapply(Wald8$p_val, Wald8$test, p.adjust, method = "bonferroni", simplify = FALSE)
  Wald8_adjusted <- Wald8
  Wald8_adjusted$p_val <- unsplit(Wald8_p_adjusted, Wald8$test)
  
  expect_equal(Wald8_adjusted, Wald9)
  
})

data(Duncan, package = "carData")
Duncan$cluster <- sample(LETTERS[1:8], size = nrow(Duncan), replace = TRUE)


test_that("Wald_test words with multiple comparisons adjustment and a single test type (Duncan data).", {
  
  Duncan_fit <- lm(prestige ~ 0 + type + income + type:income + type:education, data=Duncan)
  
  Wald5 <- Wald_test(
    Duncan_fit,
    constraints = constrain_pairwise(":education", reg_ex = TRUE),
    vcov = "CR0",
    cluster = Duncan$cluster,
    test = "HTZ"
  )
  
  Wald6 <- Wald_test(
    Duncan_fit,
    constraints = constrain_pairwise(":education", reg_ex = TRUE),
    vcov = "CR0",
    cluster = Duncan$cluster,
    test = "HTZ",
    adjustment_method = "holm"
  )
  
  Wald5_p_values <- sapply(Wald5, function(x) x$p_val) # extract p-values of Wald5
  Wald6_p_values <- sapply(Wald6, function(x) x$p_val) # extract p-values of Wald6
  
  expect_false(all(Wald5_p_values == Wald6_p_values))
  
  # get adjusted_p_values for Wald5, formatted the same as extracted p-values from Wald6
  Wald5_adjusted_p <- p.adjust(Wald5_p_values, method = "BY")
  expect_equal(Wald5_adjusted_p, Wald6_p_values)
  
  # Now using tidy = TRUE
  Wald8 <- Wald_test(
    Duncan_fit,
    constraints = constrain_pairwise(":education", reg_ex = TRUE),
    vcov = "CR0",
    cluster = Duncan$cluster,
    test = "EDF",
    adjustment_method = "none",
    tidy = TRUE
  )
  
  Wald9 <- Wald_test(
    Duncan_fit,
    constraints = constrain_pairwise(":education", reg_ex = TRUE),
    vcov = "CR0",
    cluster = Duncan$cluster,
    test = "EDF",
    adjustment_method = "fdr",
    tidy = TRUE
  )
  
  Wald8_p_adjusted <- tapply(Wald8$p_val, Wald8$test, p.adjust, method = "fdr", simplify = FALSE)
  Wald8_adjusted <- Wald8
  Wald8_adjusted$p_val <- unsplit(Wald8_p_adjusted, Wald8$test)
  
  expect_equal(Wald8_adjusted, Wald9)
  
})

test_that("Wald_test words with multiple comparisons adjustment and multiple test types  (Duncan data).", {
  
  Duncan_fit <- lm(prestige ~ 0 + type + income + type:income + type:education, data=Duncan)
    
  Wald5 <- Wald_test(
    Duncan_fit,
    constraints = constrain_pairwise(":education", reg_ex = TRUE),
    vcov = "CR0",
    cluster = Duncan$cluster,
    test = c("HTZ","chi-sq")
  )
  
  Wald6 <- Wald_test(
    Duncan_fit,
    constraints = constrain_pairwise(":education", reg_ex = TRUE),
    vcov = "CR0",
    cluster = Duncan$cluster,
    test = c("HTZ","chi-sq"),
    adjustment_method = "BY"
  )

  Wald5_p_values <- sapply(Wald5, function(x) x$p_val) # extract p-values of Wald5
  Wald6_p_values <- sapply(Wald6, function(x) x$p_val) # extract p-values of Wald6
  
  expect_false(all(Wald5_p_values == Wald6_p_values))
  
  # get adjusted_p_values for Wald5, formatted the same as extracted p-values from Wald6
  Wald5_adjusted_p <- apply(Wald5_p_values, 1, p.adjust, method = "BY", simplify = FALSE)
  Wald5_adjusted_p <- do.call(rbind, Wald5_adjusted_p)
  expect_equal(Wald5_adjusted_p, Wald6_p_values)
  
  # Now using tidy = TRUE
  
  Wald8 <- Wald_test(
    Duncan_fit,
    constraints = constrain_pairwise(":education", reg_ex = TRUE),
    vcov = "CR0",
    cluster = Duncan$cluster,
    test = c("HTA","EDF","EDT"),
    tidy = TRUE
  )
  
  Wald9 <- Wald_test(
    Duncan_fit,
    constraints = constrain_pairwise(":education", reg_ex = TRUE),
    vcov = "CR0",
    cluster = Duncan$cluster,
    test = c("HTA","EDF","EDT"),
    adjustment_method = "hommel",
    tidy = TRUE
  )
  
  Wald8_p_adjusted <- tapply(Wald8$p_val, Wald8$test, p.adjust, method = "hommel", simplify = FALSE)
  Wald8_adjusted <- Wald8
  Wald8_adjusted$p_val <- unsplit(Wald8_p_adjusted, Wald8$test)
  
  expect_equal(Wald8_adjusted, Wald9)
  
})


