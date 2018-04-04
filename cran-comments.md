## Resubmission

This is a resubmission. The new version is a maintenance release to fix compatibility issues with the car package. It also fixes the errors reported in the CRAN package check results for the previous version.

Note that the DESCRIPTION file includes a citation to Bell and McCaffrey (2002), but unfortunately that article does not have a DOI; I have included a URL instead.

## Test environments

* local Windows 7 Enterprise, R 3.4.4
* ubuntu 12.04.5 (on travis-ci), R 3.4.4, devel
* win-builder (release, devel)

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
  