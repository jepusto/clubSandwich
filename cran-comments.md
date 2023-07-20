## Resubmission

This is a re-submission in response to a notice from Kurt Hornik  on 2023-07-18. This version fixes the formatting of package version numbers in unit tests to conform to changes in `packageVersion()` in R-devel.

Please note that the DESCRIPTION file includes a citation to Bell and McCaffrey (2002), but unfortunately that article does not have a DOI; I have included a URL instead.

## Test environments

* local Windows 11 Pro, R 4.2.2
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
  linearization (8:40)
  lme (20:40)
  mlm (17:77)

  All of the identified words are spelled correctly. 

* Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1080/07350015.2016.1247004
    From: inst/doc/Wald-tests-in-clubSandwich.html
    Status: 403
    Message: Forbidden
  URL: https://doi.org/10.1257/jep.25.2.133
    From: inst/doc/panel-data-CRVE.html
    Status: 403
    Message: Forbidden
  URL: https://doi.org/10.3102/1076998615606099
    From: inst/doc/Wald-tests-in-clubSandwich.html
          inst/doc/meta-analysis-with-CRVE.html
    Status: 403
    Message: Forbidden
  URL: https://doi.org/10.4073/csr.2011.8
    From: inst/doc/meta-analysis-with-CRVE.html
    Status: 403
    Message: Forbidden
    
  The flagged URLs are correct.

Found the following (possibly) invalid DOIs:
  DOI: 10.1080/07350015.2016.1247004
    From: DESCRIPTION
    Status: Forbidden
    Message: 403
    
  The flagged DOI is correct.
  
## revdepcheck results

We checked 14 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

