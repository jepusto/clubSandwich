## Resubmission

This is a re-submission. The only change is to remove the dependency on the mathjaxr package, so that the package can be built for the Debian/GNU Linux distribution.

Please note that the DESCRIPTION file includes a citation to Bell and McCaffrey (2002), but unfortunately that article does not have a DOI; I have included a URL instead.

## Test environments

* local Windows 10 Pro, R 4.0.2
* local Windows 10 Education, R 4.0.2
* ubuntu 16.04.6 LTS (on travis-ci), R-release, R-devel
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
  URL: https://doi.org/10.1257/aer.99.4.1384
    From: man/AchievementAwardsRCT.Rd
    Status: Error
    Message: libcurl error code 56:
      	Recv failure: Connection was reset
  URL: https://doi.org/10.1257/jep.25.2.133
    From: inst/doc/panel-data-CRVE.html
    Status: Error
    Message: libcurl error code 56:
      	Recv failure: Connection was reset
  URL: https://www.jstatsoft.org/v36/i03/
    From: inst/doc/meta-analysis-with-CRVE.html
    Status: Error
    Message: libcurl error code 7:
      	Failed to connect to www.jstatsoft.org port 443: Timed out

  The flagged URLs are correct.
