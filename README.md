# r4pde

**A Companion Package for R4PDE Book**

![](https://r4pde.netlify.app/cover.png){style="float: right; padding: 20px" width="210"}I am thrilled to introduce r4pde, a comprehensive R package developed as a companion to the book "[R for Plant Disease Epidemiology](https://r4pde.netlify.app/)" (R4PDE). The companion package leverages the power of R to empower researchers and practitioners in their exploration of plant disease epidemiology. It process access to a suite of specialized R functions and data sets tailored specifically for plant disease epidemiology analysis. This package, combined with the invaluable insights and guidance offered by R4PDE, empowers you to unlock the full potential of R in tackling the complex challenges of studying plant disease epidemics.

The development version of **r4pde** is available from GitHub. The **remotes** package, available from CRAN, is required for installation.

``` r
if (!require(remotes)) {
  install.packages("remotes")
}

remotes::install_github("emdelponte/r4pde")
```

## Meta

-   Please [report any issues or bugs](https://github.com/emdelponte/r4pde/issues).

-   All code is licensed MIT

-   To cite r4pde, please use the output from citation("r4pde")
