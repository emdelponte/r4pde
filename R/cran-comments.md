## Test environments
* Local macOS (R 4.4.0)
* GitHub Actions (Ubuntu-latest, R 4.3)
* Win-builder (devel and release)

## R CMD check results

0 ERRORs | 0 WARNINGs | 2 NOTEs

* NOTE: "unable to verify current time" — this is expected due to build environment and is harmless.
* NOTE: Non-standard files/directories (‘LICENSE.md’, ‘_pkgdown.yml’, ‘docs’) — these are used for the pkgdown website and are not included in the installed package.

These notes do not require any action.

## Comments

This is the first submission of the r4pde package. It provides datasets and utility functions to support the book
*R for Plant Disease Epidemiology*, including tools for sampling design, disease progress analysis, spatial pattern evaluation,
and teaching of key epidemiological concepts.

The package depends on one Bioconductor package (`Icens`), which is listed under Imports.
