## Resubmission

This is a resubmission. 

Note that the DESCRIPTION file includes a citation to Bell and McCaffrey (2002), but unfortunately that article does not have a DOI; I have included a URL instead.

## Test environments

* local Windows 7 Enterprise, R 3.3.1
* ubuntu 12.04.5 (on travis-ci), R 3.2.5, 3.3.1, devel
* win-builder (release)

## R CMD check results

There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE

Possibly mis-spelled words in DESCRIPTION:
  Hotelling's (13:29)
  Satterthwaite (12:5)
  covariance (10:5)

  All of the identified words are spelled correctly. 

Found the following (possibly) invalid URLs:
  URL: http://dx.doi.org/10.1257/aer.99.4.1384
    From: man/AchievementAwardsRCT.Rd
    Status: Error
    Message: libcurl error code 35
    	Unknown SSL protocol error in connection to www.aeaweb.org:443
  URL: http://dx.doi.org/10.1257/jep.25.2.133
    From: inst/doc/panel-data-CRVE.html
    Status: Error
    Message: libcurl error code 35
    	Unknown SSL protocol error in connection to www.aeaweb.org:443
    	
  Both of the flagged URLs are correct.