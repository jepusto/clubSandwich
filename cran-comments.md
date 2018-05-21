## Resubmission

This is a resubmission. This submission is a maintenance release that corrects two errors identified in the "additional issues" CRAN package checks for the ATLAS and MKL BLAS. The new release adds two suggested packages (lme4 and zoo) that are used in the unit tests, but were previously not listed. The release also adds a convenience function for calculating confidence intervals based on cluster-robust standard errors and Satterthwaite degrees of freedom.

Please note that the DESCRIPTION file includes a citation to Bell and McCaffrey (2002), but unfortunately that article does not have a DOI; I have included a URL instead.

## Test environments

* local Windows 7 Enterprise, R 3.4.4
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.0, devel
* win-builder (devel, release, oldrelease)

## R CMD check results

There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE

Possibly mis-spelled words in DESCRIPTION:
  DOI (10:75)
  Hotelling's (15:49)
  McCaffrey (9:9)
  Pustejovsky (10:44)
  Satterthwaite (14:22)
  Tipton (10:60)
  eng (10:10)
  ivreg (17:21)
  mlm (16:77)
  
  All of the identified words are spelled correctly. 

Found the following (possibly) invalid URLs:
  URL: http://economics.mit.edu/faculty/angrist/data1/data/angrist
    From: man/AchievementAwardsRCT.Rd
    Status: 403
    Message: Forbidden
  URL: http://masteringmetrics.com/resources/
    From: inst/doc/panel-data-CRVE.html
    Status: 500
    Message: Internal Server Error
    
  Both of the flagged URLs are correct.
  