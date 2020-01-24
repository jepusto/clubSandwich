## Resubmission

This is a resubmission. This submission is a maintenance release that fixes unit test errors identified in the CRAN additional package checks for ATLAS and noLD.

Please note that the DESCRIPTION file includes a citation to Bell and McCaffrey (2002), but unfortunately that article does not have a DOI; I have included a URL instead.

## Test environments

* local Windows 7 Enterprise, R 3.6.0
* ubuntu 14.04.5 LTS (on travis-ci), R 3.6.2, devel
* win-builder (devel, release, oldrelease)
* r-hub:
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Ubuntu Linux 16.04 LTS, R-release, GCC
  * Fedora Linux, R-devel, clang, gfortran
  * Debian Linux, R-devel, GCC

## R CMD check results

There were no ERRORs or WARNINGs. 

There were 2 NOTES:

* Possibly mis-spelled words in DESCRIPTION:
  Hotelling's (16:49)
  McCaffrey (9:9)
  Pustejovsky (11:26)
  Satterthwaite (15:22)
  Tipton (11:42)
  lme (19:44)
  mlm (17:77)

  All of the identified words are spelled correctly. 

* Found the following (possibly) invalid URLs:
  URL: https://economics.mit.edu/faculty/angrist/data1/data/angrist
    From: man/AchievementAwardsRCT.Rd
    Status: Error
    Message: libcurl error code 60:
      	SSL certificate problem: unable to get local issuer certificate
      	(Status without verification: OK)

  The flagged URL is correct.
