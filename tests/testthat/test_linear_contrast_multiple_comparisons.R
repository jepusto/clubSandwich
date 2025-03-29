data("ChickWeight", package = "datasets")
lm_fit <- lm(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)

test_that("linear_contrast multiple comparisons p-value adjustment works correctly", {
  
  # taken from example usage
  lc1 <- linear_contrast(lm_fit,
                         vcov = "CR2",
                         cluster = ChickWeight$Chick, 
                         contrasts = constrain_pairwise("Diet.:Time", reg_ex = TRUE),
                         p_values = TRUE)
  
  lc2 <- linear_contrast(lm_fit,
                         vcov = "CR2",
                         cluster = ChickWeight$Chick,
                         contrasts = constrain_pairwise("Diet.:Time", reg_ex = TRUE),
                         p_values = TRUE, 
                         adjustment_method = "none")
  
  expect_equal(lc1, lc2) # check explicitly stating adjustment_method default doesn't change anything
  
  lc3 <- linear_contrast(lm_fit,
                         vcov = "CR2",
                         cluster = ChickWeight$Chick,
                         contrasts = constrain_pairwise("Diet.:Time", reg_ex = TRUE),
                         p_values = TRUE, 
                         adjustment_method = "holm")
  
  expect_false(all(lc2$p_val == lc3$p_val))
  expect_equal(p.adjust(lc2$p_val, method = "holm"), lc3$p_val)
  
  # test using adjustment without p_values = TRUE
  expect_warning(lc1 <- linear_contrast(lm_fit,
                                        vcov = "CR2",
                                        cluster = ChickWeight$Chick,
                                        contrasts = constrain_pairwise("Diet.:Time", reg_ex = TRUE),
                                        p_values = FALSE, adjustment_method = "hochberg"))
  lc1$p_val <- lc2$p_val
  expect_equal(lc1, lc2)
  
  # test using nonexistent adjustment type
  # changed to expect error from match.arg
  expect_error(linear_contrast(lm_fit,
                               vcov = "CR2",
                               cluster = ChickWeight$Chick, 
                               contrasts = constrain_pairwise("Diet.:Time", reg_ex = TRUE),
                               p_values = TRUE,
                               adjustment_method = "nonexistent"),
               "\'arg\' should be one of")
  
  # test using p-adjustment when results have a length of 1
  c_mat <- matrix(c(-1,1,rep(0,6)), nrow = 1)
  expect_warning(
    lc1 <- linear_contrast(lm_fit,
                           vcov = "CR2",
                           cluster = ChickWeight$Chick, 
                           contrasts = c_mat,
                           p_values = TRUE,
                           adjustment_method = "BY") # changed from BF
  )
  
  lc2 <- linear_contrast(lm_fit,
                         vcov = "CR2",
                         cluster = ChickWeight$Chick, 
                         contrasts = c_mat,
                         p_values = TRUE)
  
  expect_identical(lc1, lc2)
  
  lc4 <- linear_contrast(lm_fit,
                         vcov = "CR2",
                         cluster = ChickWeight$Chick,
                         contrasts = constrain_pairwise("Diet.$", reg_ex = TRUE),
                         p_values = TRUE, 
                         adjustment_method = "none")
  
  lc5 <- linear_contrast(lm_fit,
                         vcov = "CR2",
                         cluster = ChickWeight$Chick,
                         contrasts = constrain_pairwise("Diet.$", reg_ex = TRUE),
                         p_values = TRUE, 
                         adjustment_method = "holm")
  
  expect_false(all(lc4$p_val == lc5$p_val))
  expect_equal(p.adjust(lc4$p_val, method = "holm"), lc5$p_val)
  
  lc6 <- linear_contrast(lm_fit,
                         vcov = "CR2",
                         cluster = ChickWeight$Chick,
                         contrasts = constrain_pairwise("Diet.$", reg_ex = TRUE),
                         p_values = TRUE, 
                         adjustment_method = "hochberg")
  
  expect_false(all(lc4$p_val == lc6$p_val))
  expect_equal(p.adjust(lc4$p_val, method = "hochberg"), lc6$p_val)
  
  lc7 <- linear_contrast(lm_fit,
                         vcov = "CR2",
                         cluster = ChickWeight$Chick,
                         contrasts = constrain_pairwise("Diet.$", reg_ex = TRUE),
                         p_values = TRUE, 
                         adjustment_method = "bonferroni")
  
  expect_false(all(lc4$p_val == lc7$p_val))
  expect_equal(p.adjust(lc4$p_val, method = "bonferroni"), lc7$p_val)
  
  lc8 <- linear_contrast(lm_fit,
                         vcov = "CR2",
                         cluster = ChickWeight$Chick,
                         contrasts = constrain_pairwise("Diet.$", reg_ex = TRUE),
                         p_values = TRUE, 
                         adjustment_method = "BH")
  
  expect_false(all(lc4$p_val == lc8$p_val))
  expect_equal(p.adjust(lc4$p_val, method = "BH"), lc8$p_val)
  
  lc9 <- linear_contrast(lm_fit,
                         vcov = "CR2",
                         cluster = ChickWeight$Chick,
                         contrasts = constrain_pairwise("Diet.$", reg_ex = TRUE),
                         p_values = TRUE, 
                         adjustment_method = "BY")
  
  expect_false(all(lc4$p_val == lc9$p_val))
  expect_equal(p.adjust(lc4$p_val, method = "BY"), lc9$p_val)
  
})
