#' Achievement Awards Demonstration program
#' 
#' Data from a randomized trial of the Achievement Awards
#' Demonstration program, reported in Angrist & Lavy (2009).
#' 
#' @format A data frame with 16526 rows and 21 variables: \describe{ 
#'   \item{school_id}{Fictitious school identification number}
#'   \item{school_type}{Factor identifying the school type (Arab religious, Jewish religious, Jewish secular)}
#'   \item{pair}{Number of treatment pair. Note that 7 is a triple.} 
#'   \item{treated}{Indicator for whether school was in treatment group}
#'   \item{year}{Cohort year}
#'   \item{student_id}{Fictitious student identification number}
#'   \item{sex}{Factor identifying student sex}
#'   \item{siblings}{Number of siblings}
#'   \item{immigrant}{Indicator for immigrant status}
#'   \item{father_ed}{Father's level of education}
#'   \item{mother_ed}{Mother's level of education}
#'   \item{Bagrut_status}{Indicator for Bagrut attainment}
#'   \item{attempted}{Number of Bagrut units attempted}
#'   \item{awarded}{Number of Bagrut units awarded}
#'   \item{achv_math}{Indicator for satisfaction of math requirement}
#'   \item{achv_english}{Indicator for satisfaction of English requirement}
#'   \item{achv_hebrew}{Indicator for satisfaction of Hebrew requirement}
#'   \item{lagscore}{Lagged Bagrut score}
#'   \item{qrtl}{Quartile within distribution of lagscore, calculated by cohort and sex}
#'   \item{half}{Lower or upper half within distribution of lagscore, calculated by cohort and sex}
#'   }
#'   
#' @source \href{https://economics.mit.edu/people/faculty/josh-angrist/angrist-data-archive}{Angrist Data Archive}
#'   
#' @references Angrist, J. D., & Lavy, V. (2009). The effects of high stakes 
#'   high school achievement awards : Evidence from a randomized trial.
#'   \emph{American Economic Review, 99}(4), 1384-1414.
#'   \doi{10.1257/aer.99.4.1384}
#'   

"AchievementAwardsRCT"



#' Dropout prevention/intervention program effects
#'
#' A dataset containing estimated effect sizes, variances, and covariates from a
#' meta-analysis of dropout prevention/intervention program effects, conducted
#' by Wilson et al. (2011). Missing observations were imputed.
#'
#' @format A data frame with 385 rows and 18 variables: \describe{
#'   \item{LOR1}{log-odds ratio measuring the intervention effect}
#'   \item{varLOR}{estimated sampling variance of the log-odds ratio}
#'   \item{studyID}{unique identifier for each study} \item{studySample}{unique
#'   identifier for each sample within a study} \item{study_design}{study design
#'   (randomized, matched, or non-randomized and unmatched)}
#'   \item{outcome}{outcome measure for the intervention effect is estimated
#'   (school dropout, school enrollment, graduation, graduation or GED receipt)}
#'   \item{evaluator_independence}{degree of evaluator independence
#'   (independent, indirect but influential, involved in planning but not
#'   delivery, involved in delivery)} \item{implementation_quality}{level of
#'   implementation quality (clear problems, possible problems, no apparent
#'   problems)} \item{program_site}{Program delivery site (community, mixed,
#'   school classroom, school but outside of classroom)}
#'   \item{attrition}{Overall attrition (proportion)}
#'   \item{group_equivalence}{pretest group-equivalence log-odds ratio}
#'   \item{adjusted}{adjusted or unadjusted data used to calculate intervention
#'   effect} \item{male_pct}{proportion of the sample that is male}
#'   \item{white_pct}{proportion of the sample that is white}
#'   \item{average_age}{average age of the sample} \item{duration}{program
#'   duration (in weeks)} \item{service_hrs}{program contact hours per week}
#'   \item{big_study}{indicator for the 32 studies with 3 or more effect sizes}
#'   }
#'
#' @source Wilson, S. J., Lipsey, M. W., Tanner-Smith, E., Huang, C. H., &
#'   Steinka-Fry, K. T. (2011). Dropout prevention and intervention programs:
#'   Effects on school completion and dropout Among school-aged children and
#'   youth: A systematic review. _Campbell Systematic Reviews, 7_(1), 1-61.
#'   \doi{10.4073/csr.2011.8}
#'
#' @references Wilson, S. J., Lipsey, M. W., Tanner-Smith, E., Huang, C. H., &
#'   Steinka-Fry, K. T. (2011). Dropout prevention and intervention programs:
#'   Effects on school completion and dropout Among school-aged children and
#'   youth: A systematic review. _Campbell Systematic Reviews, 7_(1), 1-61.
#'   \doi{10.4073/csr.2011.8}
#'
#'   Tipton, E., & Pustejovsky, J. E. (2015). Small-sample adjustments for tests
#'   of moderators and model fit using robust variance estimation in
#'   meta-regression. _Journal of Educational and Behavioral Statistics, 40_(6), 604-634.
#'   \doi{10.3102/1076998615606099}
#'   

"dropoutPrevention"


#' State-level annual mortality rates by cause among 18-20 year-olds
#' 
#' A dataset containing state-level annual mortality rates for select causes of
#' death, as well as data related to the minimum legal drinking age and alcohol
#' consumption.
#' 
#' @format A data frame with 5508 rows and 12 variables: \describe{ 
#'   \item{year}{Year of observation} 
#'   \item{state}{identifier for state} 
#'   \item{count}{Number of deaths} 
#'   \item{pop}{Population size} 
#'   \item{legal}{Proportion of 18-20 year-old population that is legally allowed to drink} 
#'   \item{beertaxa}{Beer taxation rate} 
#'   \item{beerpercap}{Beer consumption per capita} 
#'   \item{winepercap}{Wine consumption per capita} 
#'   \item{spiritpercap}{Spirits consumption per capita} 
#'   \item{totpercap}{Total alcohol consumption per capita} 
#'   \item{mrate}{Mortality rate per 10,000} 
#'   \item{cause}{Cause of death} 
#'   }
#'   
#' @source
#'   \href{http://masteringmetrics.com/wp-content/uploads/2015/01/deaths.dta}{Mastering
#'   'Metrics data archive}
#'   
#' @references
#' 
#' Angrist, J. D., and Pischke, J. S. (2014). _Mastering'metrics: the path from
#' cause to effect_. Princeton University Press, 2014.
#' 
#' Carpenter, C., & Dobkin, C. (2011). The minimum legal drinking age and public
#' health. _Journal of Economic Perspectives, 25_(2), 133-156.
#' \doi{10.1257/jep.25.2.133}
#' 
 
"MortalityRates"

#' Randomized experiments on SAT coaching
#' 
#' Effect sizes from studies on the effects of SAT coaching,
#' reported in Kalaian and Raudenbush (1996)
#' 
#' @format A data frame with 67 rows and 11 variables: 
#' \describe{ 
#'   \item{study}{Study identifier}
#'   \item{year}{Year of publication} 
#'   \item{test}{Character string indicating whether effect size corresponds to outcome on verbal (SATV) or math (SATM) test}
#'   \item{d}{Effect size estimate (Standardized mean difference)} 
#'   \item{V}{Variance of effect size estimate} 
#'   \item{nT}{Sample size in treatment condition} 
#'   \item{nC}{Sample size in control condition} 
#'   \item{study_type}{Character string indicating whether study design used a matched, non-equivalent, or randomized control group} 
#'   \item{hrs}{Hours of coaching} 
#'   \item{ETS}{Indicator variable for Educational Testing Service} 
#'   \item{homework}{Indicator variable for homework} 
#'   }
#'   
#' @references Kalaian, H. A. & Raudenbush, S. W. (1996). A multivariate mixed 
#'   linear model for meta-analysis. \emph{Psychological Methods, 1}(3),
#'   227-235. 
#'   \doi{10.1037/1082-989X.1.3.227}
#'   

"SATcoaching"
