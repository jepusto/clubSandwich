## Resubmission

This is a re-submission. This version includes miscellaneous bug fixes, error-handling, and corrections to the documentation.

Please note that the DESCRIPTION file includes a citation to Bell and McCaffrey (2002), but unfortunately that article does not have a DOI; I have included a URL instead.

## Test environments

* local Windows 11 Pro, R 4.1.2
* ubuntu 20.04.3 LTS (on Github), R devel, release, oldrelease
* macOS-latest (on Github), R release
* windows-latest (on Github), R release
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
  URL: https://doi.org/10.3102/1076998615606099
    From: man/dropoutPrevention.Rd
    Status: 503
    Message: Service Unavailable
  URL: https://doi.org/10.1037/1082-989X.1.3.227
    From: man/SATcoaching.Rd
    Status: 400
    Message: Bad Request
  URL: https://doi.org/10.1037/met0000011
    From: inst/doc/meta-analysis-with-CRVE.html
    Status: 400
    Message: Bad Request
  URL: https://economics.mit.edu/faculty/angrist/data1/data/angrist
    From: man/AchievementAwardsRCT.Rd
    Status: Error
    Message: SSL certificate problem: unable to get local issuer certificate
    
  The flagged URLs are correct.

## revdepcheck results

We checked 7 reverse dependencies (6 from CRAN + 1 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 1 packages

Issues with CRAN packages are summarised below.

### Failed to check

* ggeffects (NA)
