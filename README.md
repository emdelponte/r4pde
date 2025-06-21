# r4pde

**r4pde** is an R package developed as a companion to the book [*R for Plant Disease Epidemiology*](https://r4pde.netlify.app/).  
It provides access to a suite of specialized R functions and datasets tailored for teaching and analyzing plant disease epidemiology.

This package supports researchers, students, and practitioners by offering tools for sampling design, disease progress analysis, spatial epidemiology, and related applications.

## Installation

The development version of **r4pde** is available from GitHub. To install it along with its dependencies (including Bioconductor packages), use the [`pak`](https://pak.r-lib.org/) package:

```r
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pkg_install("emdelponte/r4pde")
```

Alternatively, using **remotes**:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("emdelponte/r4pde")
```

> Note: The package depends on `Icens`, which is available from Bioconductor. Both `pak` and `remotes` will install it automatically if needed.

## Reporting Issues

Please [report any bugs or issues](https://github.com/emdelponte/r4pde/issues) via the GitHub issue tracker.

## License and Citation

- This package is released under the MIT license.
- To cite **r4pde**, please use:

```r
citation("r4pde")
```
