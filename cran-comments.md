## Resubmission

This is a resubmission. This submission is a maintenance release that corrects an error identified in the CRAN package checks for debian-gcc. The release also updates some internals for consistency with changes in the metafor package.

Please note that the DESCRIPTION file includes a citation to Bell and McCaffrey (2002), but unfortunately that article does not have a DOI; I have included a URL instead.

## Test environments

* local Windows 7 Enterprise, R 3.6.0
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.2, devel
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
  