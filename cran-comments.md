## Resubmission

This is a resubmission. In this version I have:

* Added epiR to Suggests (used in tests).
* Made the sad-validation vignette robust to non-convergence of an optional model refit.
* Moved terra from Imports to Suggests, since it is not available on the r-devel Debian check machines (current 0.2.2 ERROR).

## R CMD check results

0 errors | 0 warnings | 2 notes
