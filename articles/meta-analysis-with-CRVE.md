# Meta-analysis with cluster-robust variance estimation

This vignette demonstrates how to use the `clubSandwich` package to
conduct a meta-analysis of dependent effect sizes with robust variance
estimation. Tests of meta-regression coefficients and F-tests of
multiple-coefficient hypotheses are calculated using small-sample
corrections proposed by Tipton (2015) and Tipton and Pustejovsky (2015).
The example uses a dataset of effect sizes from a Campbell Collaboration
systematic review of dropout prevention programs, conducted by Sandra Jo
Wilson and colleagues (2011).

The original analysis included a meta-regression with covariates that
capture methodological, participant, and program characteristics. The
regression specification used here is similar to Model III from Wilson
et al. (2011), but treats the `evaluator_independence` and
`implementation_quality` variables as categorical rather than
interval-level. Also, the original analysis clustered at the level of
the sample (some studies reported results from multiple samples),
whereas here we cluster at the study level. The meta-regression can be
fit in several different ways. We first demonstrate using the `robumeta`
package (Fisher & Tipton, 2015) and then using the `metafor` package
(Viechtbauer, 2010).

## robumeta model

``` r

library(clubSandwich)
library(robumeta)
data(dropoutPrevention)

# clean formatting
names(dropoutPrevention)[7:8] <- c("eval","implement")
levels(dropoutPrevention$eval) <- c("independent","indirect","planning","delivery")
levels(dropoutPrevention$implement) <- c("low","medium","high")
levels(dropoutPrevention$program_site) <- c("community","mixed","classroom","school")
levels(dropoutPrevention$study_design) <- c("matched","unmatched","RCT")
levels(dropoutPrevention$adjusted) <- c("no","yes")

m3_robu <- robu(LOR1 ~ study_design + attrition + group_equivalence + adjusted
                + outcome + eval + male_pct + white_pct + average_age
                + implement + program_site + duration + service_hrs, 
                data = dropoutPrevention, studynum = studyID, var.eff.size = varLOR, 
                modelweights = "HIER")
print(m3_robu)
```

    ## RVE: Hierarchical Effects Model with Small-Sample Corrections 
    ## 
    ## Model: LOR1 ~ study_design + attrition + group_equivalence + adjusted + outcome + eval + male_pct + white_pct + average_age + implement + program_site + duration + service_hrs 
    ## 
    ## Number of clusters = 152 
    ## Number of outcomes = 385 (min = 1 , mean = 2.53 , median = 1 , max = 30 )
    ## Omega.sq = 0.24907 
    ## Tau.sq = 0.1024663 
    ## 
    ##                           Estimate   StdErr t-value  dfs    P(|t|>) 95% CI.L 95% CI.U Sig
    ## 1           X.Intercept.  0.016899 0.615399  0.0275 16.9 0.97841541 -1.28228  1.31608    
    ## 2  study_designunmatched -0.002626 0.185142 -0.0142 40.5 0.98875129 -0.37667  0.37141    
    ## 3        study_designRCT -0.086872 0.140044 -0.6203 38.6 0.53869676 -0.37024  0.19650    
    ## 4              attrition  0.118889 0.247228  0.4809 15.5 0.63732597 -0.40666  0.64444    
    ## 5      group_equivalence  0.502463 0.195838  2.5657 28.7 0.01579282  0.10174  0.90318  **
    ## 6            adjustedyes -0.322480 0.125413 -2.5713 33.8 0.01470796 -0.57741 -0.06755  **
    ## 7        outcomeenrolled  0.097059 0.139842  0.6941 16.5 0.49727848 -0.19862  0.39274    
    ## 8      outcomegraduation  0.147643 0.134938  1.0942 30.2 0.28253825 -0.12786  0.42315    
    ## 9  outcomegraduation.ged  0.258034 0.169134  1.5256 16.3 0.14632629 -0.10006  0.61613    
    ## 10          evalindirect -0.765085 0.399109 -1.9170  6.2 0.10212896 -1.73406  0.20389    
    ## 11          evalplanning -0.920874 0.346536 -2.6574  5.6 0.04027061 -1.78381 -0.05794  **
    ## 12          evaldelivery -0.916673 0.304303 -3.0124  4.7 0.03212299 -1.71432 -0.11903  **
    ## 13              male_pct  0.167965 0.181538  0.9252 16.4 0.36824526 -0.21609  0.55202    
    ## 14             white_pct  0.022915 0.149394  0.1534 21.8 0.87950385 -0.28704  0.33287    
    ## 15           average_age  0.037102 0.027053  1.3715 21.2 0.18458247 -0.01913  0.09333    
    ## 16       implementmedium  0.411779 0.128898  3.1946 26.7 0.00358205  0.14714  0.67642 ***
    ## 17         implementhigh  0.658570 0.123874  5.3164 34.6 0.00000635  0.40699  0.91015 ***
    ## 18     program_sitemixed  0.444384 0.172635  2.5741 28.6 0.01550504  0.09109  0.79768  **
    ## 19 program_siteclassroom  0.426658 0.159773  2.6704 37.4 0.01115192  0.10303  0.75028  **
    ## 20    program_siteschool  0.262517 0.160519  1.6354 30.1 0.11236814 -0.06525  0.59028    
    ## 21              duration  0.000427 0.000873  0.4895 36.7 0.62736846 -0.00134  0.00220    
    ## 22           service_hrs -0.003434 0.005012 -0.6852 36.7 0.49752503 -0.01359  0.00672    
    ## ---
    ## Signif. codes: < .01 *** < .05 ** < .10 *
    ## ---
    ## Note: If df < 4, do not trust the results

Note that `robumeta` produces small-sample corrected standard errors and
t-tests, and so there is no need to repeat those calculations with
`clubSandwich`. The `eval` variable has four levels, and it might be of
interest to test whether the average program effects differ by the
degree of evaluator independence. The null hypothesis in this case is
that the 10th, 11th, and 12th regression coefficients are all equal to
zero. A small-sample adjusted F-test for this hypothesis can be obtained
as follows. The `vcov = "CR2"` option means that the standard errors
will be corrected using the bias-reduced linearization estimator
described in Tipton and Pustejovsky (2015).

``` r

Wald_test(m3_robu, constraints = constrain_zero(10:12), vcov = "CR2")
```

    ##  test Fstat df_num df_denom  p_val sig
    ##   HTZ  2.78      3     16.8 0.0732   .

By default, the `Wald_test` function provides an F-type test with
degrees of freedom estimated using the approximate Hotelling’s $`T^2_Z`$
method. The test has less than 17 degrees of freedom, even though there
are 152 independent studies in the data, and has a p-value that is not
quite significant at conventional levels. The low degrees of freedom are
a consequence of the fact that one of the levels of
`evaluator independence` has only a few effect sizes in it:

``` r

table(dropoutPrevention$eval)
```

    ## 
    ## independent    indirect    planning    delivery 
    ##           6          33          43         303

## metafor model

`clubSandwich` also works with models fit using the `metafor` package.
Here we re-fit the same regression specification, but use REML to
estimate the variance components (`robumeta` uses a method-of-moments
estimator), as well as a somewhat different weighting scheme than that
used in `robumeta`.

``` r

library(metafor)
m3_metafor <- rma.mv(LOR1 ~ study_design + attrition + group_equivalence + adjusted
                      + outcome + eval
                      + male_pct + white_pct + average_age
                      + implement + program_site + duration + service_hrs, 
                      V = varLOR, random = list(~ 1 | studyID, ~ 1 | studySample),
                     data = dropoutPrevention)
summary(m3_metafor)
```

    ## 
    ## Multivariate Meta-Analysis Model (k = 385; method: REML)
    ## 
    ##    logLik   Deviance        AIC        BIC       AICc   
    ## -489.0357   978.0714  1026.0714  1119.5371  1029.6217   
    ## 
    ## Variance Components:
    ## 
    ##             estim    sqrt  nlvls  fixed       factor 
    ## sigma^2.1  0.2274  0.4769    152     no      studyID 
    ## sigma^2.2  0.1145  0.3384    317     no  studySample 
    ## 
    ## Test for Residual Heterogeneity:
    ## QE(df = 363) = 1588.4397, p-val < .0001
    ## 
    ## Test of Moderators (coefficients 2:22):
    ## QM(df = 21) = 293.8694, p-val < .0001
    ## 
    ## Model Results:
    ## 
    ##                        estimate      se     zval    pval    ci.lb    ci.ub      
    ## intrcpt                  0.5296  0.7250   0.7304  0.4651  -0.8915   1.9506      
    ## study_designunmatched   -0.0494  0.1722  -0.2871  0.7741  -0.3870   0.2881      
    ## study_designRCT          0.0653  0.1628   0.4010  0.6884  -0.2538   0.3843      
    ## attrition               -0.1366  0.2429  -0.5623  0.5739  -0.6126   0.3395      
    ## group_equivalence        0.4071  0.1573   2.5877  0.0097   0.0988   0.7155   ** 
    ## adjustedyes             -0.3581  0.1532  -2.3371  0.0194  -0.6585  -0.0578    * 
    ## outcomeenrolled         -0.2831  0.0771  -3.6709  0.0002  -0.4343  -0.1320  *** 
    ## outcomegraduation       -0.0913  0.0657  -1.3896  0.1646  -0.2201   0.0375      
    ## outcomegraduation/ged    0.6983  0.0805   8.6750  <.0001   0.5406   0.8561  *** 
    ## evalindirect            -0.7530  0.4949  -1.5214  0.1282  -1.7230   0.2171      
    ## evalplanning            -0.7700  0.4869  -1.5814  0.1138  -1.7242   0.1843      
    ## evaldelivery            -1.0016  0.4600  -2.1774  0.0294  -1.9033  -0.1000    * 
    ## male_pct                 0.1021  0.1715   0.5951  0.5518  -0.2341   0.4382      
    ## white_pct                0.1223  0.1804   0.6777  0.4979  -0.2313   0.4758      
    ## average_age              0.0061  0.0291   0.2091  0.8344  -0.0509   0.0631      
    ## implementmedium          0.4738  0.1609   2.9445  0.0032   0.1584   0.7892   ** 
    ## implementhigh            0.6318  0.1471   4.2965  <.0001   0.3436   0.9201  *** 
    ## program_sitemixed        0.3289  0.2413   1.3631  0.1729  -0.1440   0.8019      
    ## program_siteclassroom    0.2920  0.1736   1.6821  0.0926  -0.0482   0.6321    . 
    ## program_siteschool       0.1616  0.1898   0.8515  0.3945  -0.2104   0.5337      
    ## duration                 0.0013  0.0009   1.3423  0.1795  -0.0006   0.0031      
    ## service_hrs             -0.0003  0.0047  -0.0654  0.9478  -0.0096   0.0090      
    ## 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

`metafor` produces model-based standard errors, t-tests, and confidence
intervals. The `coef_test` function from `clubSandwich` will calculate
robust standard errors and robust t-tests for each of the coefficients:

``` r

coef_test(m3_metafor, vcov = "CR2")
```

    ## Alternative hypothesis: two-sided 
    ##                  Coef.  Estimate       SE Null value  t-stat d.f. (Satt) p-val (Satt) Sig.
    ##                intrcpt  0.529569 0.724851          0  0.7306       20.08      0.47347     
    ##  study_designunmatched -0.049434 0.204152          0 -0.2421       58.42      0.80952     
    ##        study_designRCT  0.065272 0.149146          0  0.4376       53.17      0.66342     
    ##              attrition -0.136575 0.306429          0 -0.4457       10.52      0.66485     
    ##      group_equivalence  0.407108 0.210917          0  1.9302       23.10      0.06595    .
    ##            adjustedyes -0.358124 0.136132          0 -2.6307       43.20      0.01176    *
    ##        outcomeenrolled -0.283124 0.237199          0 -1.1936        7.08      0.27108     
    ##      outcomegraduation -0.091295 0.091465          0 -0.9981        9.95      0.34188     
    ##  outcomegraduation/ged  0.698328 0.364882          0  1.9138        8.02      0.09188    .
    ##           evalindirect -0.752994 0.447670          0 -1.6820        6.56      0.13929     
    ##           evalplanning -0.769968 0.403898          0 -1.9063        6.10      0.10446     
    ##           evaldelivery -1.001648 0.355989          0 -2.8137        4.89      0.03834    *
    ##               male_pct  0.102055 0.148410          0  0.6877        9.68      0.50782     
    ##              white_pct  0.122255 0.141470          0  0.8642       16.88      0.39961     
    ##            average_age  0.006084 0.033387          0  0.1822       15.79      0.85772     
    ##        implementmedium  0.473789 0.148660          0  3.1871       22.44      0.00419   **
    ##          implementhigh  0.631842 0.138073          0  4.5761       28.68      < 0.001  ***
    ##      program_sitemixed  0.328941 0.196848          0  1.6710       27.47      0.10607     
    ##  program_siteclassroom  0.291952 0.146014          0  1.9995       42.70      0.05195    .
    ##     program_siteschool  0.161640 0.171700          0  0.9414       29.27      0.35420     
    ##               duration  0.001270 0.000978          0  1.2988       31.96      0.20332     
    ##            service_hrs -0.000309 0.004828          0 -0.0641       49.63      0.94915

Note that `coef_test` assumed that it should cluster based on `studyID`,
which is the outer-most random effect in the metafor model. This can be
specified explicitly by including the option
`cluster = dropoutPrevention$studyID` in the call.

The F-test for degree of evaluator independence uses the same syntax as
before:

``` r

Wald_test(m3_metafor, constraints = constrain_zero(10:12), vcov = "CR2")
```

    ##  test Fstat df_num df_denom  p_val sig
    ##   HTZ  2.71      3     18.3 0.0753   .

Despite some differences in weighting schemes, the p-value is very close
to the result obtained using `robumeta`.

## References

Fisher, Z., & Tipton, E. (2015). robumeta: An R-package for robust
variance estimation in meta-analysis.
[arXiv:1503.02220](https://arxiv.org/abs/1503.02220)

Tipton, E. (2015). Small sample adjustments for robust variance
estimation with meta-regression. *Psychological Methods, 20*(3),
375-393. <https://doi.org/10.1037/met0000011>

Tipton, E., & Pustejovsky, J. E. (2015). Small-sample adjustments for
tests of moderators and model fit using robust variance estimation in
meta-regression. *Journal of Educational and Behavioral Statistics,
40*(6), 604-634. <https://doi.org/10.3102/1076998615606099>

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor
package. *Journal of Statistical Software, 36*(3), 1-48. URL:
<https://doi.org/10.18637/jss.v036.i03>

Wilson, S. J., Lipsey, M. W., Tanner-Smith, E., Huang, C. H., &
Steinka-Fry, K. T. (2011). Dropout prevention and intervention programs:
Effects on school completion and dropout Among school-aged children and
youth: A systematic review. *Campbell Systematic Reviews, 7*(1), 1-61.
<https://doi.org/10.4073/csr.2011.8>
