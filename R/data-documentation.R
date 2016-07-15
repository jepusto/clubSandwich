#' Dropout prevention/intervention program effects
#' 
#' A dataset containing estimated effect sizes, variances, and covariates from a
#' meta-analysis of dropout prevention/intervention program effects, conducted
#' by Wilson et al. (2011). Missing observations were imputed. 
#' 
#' @format A data frame with 385 rows and 18 variables: \describe{ 
#'   \item{LOR1}{log-odds ratio measuring the intervention effect} 
#'   \item{varLOR}{estimated sampling variance of the log-odds ratio} 
#'   \item{studyID}{unique identifier for each study}
#'   \item{studySample}{unique identifier for each sample within a study}
#'   \item{study_design}{study design (randomized, matched, or non-randomized and unmatched)}
#'   \item{outcome}{outcome measure for the intervention effect is estimated (school dropout, 
#'   school enrollment, graduation, graduation or GED receipt)}
#'   \item{evaluator_independence}{degree of evaluator independence (independent, indirect 
#'   but influential, involved in planning but not delivery, involved in delivery)}
#'   \item{implementation_quality}{level of implementation quality (clear problems, 
#'   possible problems, no apparent problems)}
#'   \item{program_site}{Program delivery site (community, mixed, school classroom, 
#'   school but outside of classroom)}
#'   \item{attrition}{Overall attrition (proportion)}
#'   \item{group_equivalence}{pretest group-equivalence log-odds ratio}
#'   \item{adjusted}{adjusted or unadjusted data used to calculate intervention effect}
#'   \item{male_pct}{proportion of the sample that is male}
#'   \item{white_pct}{proportion of the sample that is white}
#'   \item{average_age}{average age of the sample}
#'   \item{duration}{program duration (in weeks)}
#'   \item{service_hrs}{program contact hours per week}
#'   \item{big_study}{indicator for the 32 studies with 3 or more effect sizes}
#'   }
#'   
#' @source Wilson, S. J., Lipsey, M. W., Tanner-Smith, E., Huang, C. H., &
#'   Steinka-Fry, K. T. (2011). Dropout prevention and intervention programs:
#'   Effects on school completion and dropout Among school-aged children and
#'   youth: A systematic review. Campbell Systematic Reviews, 7(8).
#'   
#' @references Wilson, S. J., Lipsey, M. W., Tanner-Smith, E., Huang, C. H., &
#'   Steinka-Fry, K. T. (2011). Dropout prevention and intervention programs:
#'   Effects on school completion and dropout Among school-aged children and
#'   youth: A systematic review. Campbell Systematic Reviews, 7(8).
#'   
#'   Tipton, E., & Pustejovsky, J. E. (2015). Small-sample adjustments for tests
#'   of moderators and model fit using robust variance estimation in
#'   meta-regression.
#'   

"dropoutPrevention"


#' Achievement Awards Demonstration program
#' 
#' A dataset containing data from a randomized trial of the Achievement Awards
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
#' @source \href{http://economics.mit.edu/faculty/angrist/data1/data/angrist}{Angrist Data Archive}
#'   
#' @references Angrist, J. D., & Lavy, V. (2009). The effects of high stakes 
#'   high school achievement awards : Evidence from a randomized trial.
#'   \emph{American Economic Review, 99}(4), 1384-1414.
#'   doi:\href{http://dx.doi.org/10.1257/aer.99.4.1384}{10.1257/aer.99.4.1384}
#'   

"AchievementAwardsRCT"
