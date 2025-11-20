# Randomized experiments on SAT coaching

Effect sizes from studies on the effects of SAT coaching, reported in
Kalaian and Raudenbush (1996)

## Usage

``` r
SATcoaching
```

## Format

A data frame with 67 rows and 11 variables:

- study:

  Study identifier

- year:

  Year of publication

- test:

  Character string indicating whether effect size corresponds to outcome
  on verbal (SATV) or math (SATM) test

- d:

  Effect size estimate (Standardized mean difference)

- V:

  Variance of effect size estimate

- nT:

  Sample size in treatment condition

- nC:

  Sample size in control condition

- study_type:

  Character string indicating whether study design used a matched,
  non-equivalent, or randomized control group

- hrs:

  Hours of coaching

- ETS:

  Indicator variable for Educational Testing Service

- homework:

  Indicator variable for homework

## References

Kalaian, H. A. & Raudenbush, S. W. (1996). A multivariate mixed linear
model for meta-analysis. *Psychological Methods, 1*(3), 227-235.
[doi:10.1037/1082-989X.1.3.227](https://doi.org/10.1037/1082-989X.1.3.227)
