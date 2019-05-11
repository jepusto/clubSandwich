## Resubmission

This is a resubmission. This submission is a maintenance release that corrects errors identified in the CRAN package checks for debian-gcc and the additional checks for ATLAS, MKL, and noLD. The release also includes a bug fix for methods supporting the plm package and some minor tweaks to displayed output.

Please note that the DESCRIPTION file includes a citation to Bell and McCaffrey (2002), but unfortunately that article does not have a DOI; I have included a URL instead.

## Test environments

* local Windows 7 Enterprise, R 3.6.0
* ubuntu 14.04.5 LTS (on travis-ci), R 3.6.0, devel
* win-builder (devel, release, oldrelease)
* r-hub:
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Ubuntu Linux 16.04 LTS, R-release, GCC
  * Fedora Linux, R-devel, clang, gfortran
  * Debian Linux, R-devel, GCC

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
  URL: https://economics.mit.edu/faculty/angrist/data1/data/angrist
    From: man/AchievementAwardsRCT.Rd
    Status: Error
    Message: libcurl error code 60:
      	SSL certificate problem: unable to get local issuer certificate
      	(Status without verification: OK)

  The flagged URL is correct.
  