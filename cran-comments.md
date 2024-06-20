## Resubmission

This is a re-submission. This version is a maintenance release. The main change is to deprecate two helper functions that have been made obsolete by updates to other packages. We have also corrected formatting issues in the .Rd files that were noted in CRAN checks of the previous release.

Please note that the DESCRIPTION file includes a citation to Bell and McCaffrey (2002), but unfortunately that article does not have a DOI; I have included a URL instead.

## Test environments

* local Windows 11 Pro, R 4.3.3
* ubuntu 20.04.3 LTS (on Github), R devel, release, oldrelease
* macOS-latest (on Github), R release
* windows-latest (on Github), R release
* win-builder (devel, release, oldrelease)
* mac-builder (release)

## R CMD check results

There were no ERRORs, WARNINGs, or NOTEs. 

## recheck results

We checked 15 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

