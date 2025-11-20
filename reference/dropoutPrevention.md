# Dropout prevention/intervention program effects

A dataset containing estimated effect sizes, variances, and covariates
from a meta-analysis of dropout prevention/intervention program effects,
conducted by Wilson et al. (2011). Missing observations were imputed.

## Usage

``` r
dropoutPrevention
```

## Format

A data frame with 385 rows and 18 variables:

- LOR1:

  log-odds ratio measuring the intervention effect

- varLOR:

  estimated sampling variance of the log-odds ratio

- studyID:

  unique identifier for each study

- studySample:

  unique identifier for each sample within a study

- study_design:

  study design (randomized, matched, or non-randomized and unmatched)

- outcome:

  outcome measure for the intervention effect is estimated (school
  dropout, school enrollment, graduation, graduation or GED receipt)

- evaluator_independence:

  degree of evaluator independence (independent, indirect but
  influential, involved in planning but not delivery, involved in
  delivery)

- implementation_quality:

  level of implementation quality (clear problems, possible problems, no
  apparent problems)

- program_site:

  Program delivery site (community, mixed, school classroom, school but
  outside of classroom)

- attrition:

  Overall attrition (proportion)

- group_equivalence:

  pretest group-equivalence log-odds ratio

- adjusted:

  adjusted or unadjusted data used to calculate intervention effect

- male_pct:

  proportion of the sample that is male

- white_pct:

  proportion of the sample that is white

- average_age:

  average age of the sample

- duration:

  program duration (in weeks)

- service_hrs:

  program contact hours per week

- big_study:

  indicator for the 32 studies with 3 or more effect sizes

## Source

Wilson, S. J., Lipsey, M. W., Tanner-Smith, E., Huang, C. H., &
Steinka-Fry, K. T. (2011). Dropout prevention and intervention programs:
Effects on school completion and dropout Among school-aged children and
youth: A systematic review. \_Campbell Systematic Reviews, 7\_(1), 1-61.
[doi:10.4073/csr.2011.8](https://doi.org/10.4073/csr.2011.8)

## References

Wilson, S. J., Lipsey, M. W., Tanner-Smith, E., Huang, C. H., &
Steinka-Fry, K. T. (2011). Dropout prevention and intervention programs:
Effects on school completion and dropout Among school-aged children and
youth: A systematic review. \_Campbell Systematic Reviews, 7\_(1), 1-61.
[doi:10.4073/csr.2011.8](https://doi.org/10.4073/csr.2011.8)

Tipton, E., & Pustejovsky, J. E. (2015). Small-sample adjustments for
tests of moderators and model fit using robust variance estimation in
meta-regression. \_Journal of Educational and Behavioral Statistics,
40\_(6), 604-634.
[doi:10.3102/1076998615606099](https://doi.org/10.3102/1076998615606099)
