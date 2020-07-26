## Resubmission

This is a resubmission. This submission adds new functionality for Wald tests, along with expanded documentation and a vignette. This submission also fixes unit test errors that occur with MKL and noLD builds. 

Please note that the DESCRIPTION file includes a citation to Bell and McCaffrey (2002), but unfortunately that article does not have a DOI; I have included a URL instead.

## Test environments

* local Windows 7 Enterprise, R 4.0.2
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

  The flagged URLs are correct.
