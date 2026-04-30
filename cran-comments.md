## Resubmission

This release adds support for two more classes of model objects, including iv_robust() models from the estimatr package and mmrm() models from the mmrm package. 
It also fixes a bug in an internal function introduced by changes to error message structure in R-devel, which corrects all CRAN checks that are currently failing.

Please note that the DESCRIPTION file includes a citation to Bell and McCaffrey (2002), but unfortunately that article does not have a DOI; I have included a URL instead.

## Test environments

* local Windows 11 Pro, R 4.4.3
* ubuntu 20.04.3 LTS (on Github), R devel, release, oldrelease
* macOS-latest (on Github), R release
* windows-latest (on Github), R release
* win-builder (devel, release, oldrelease)
* mac-builder (release)

## R CMD check results

There were no ERRORs, no WARNINGs, and no NOTEs.

## recheck results

We checked 21 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

