# Plot the output of temporal compartmental epidemic models

Helper used in the chapters on temporal compartmental models of the book
*R for Plant Disease Epidemiology* (R4PDE). It draws the dynamics of the
host compartments simulated by HLIR-type models, and the disease
progress and inoculum level across seasons of a multi-seasonal epidemic.

## Usage

``` r
plot_model_results(
  data,
  type = c("HLIR", "HLIRP", "HLIRP (with host growth)", "MULTI"),
  print = TRUE
)
```

## Arguments

- data:

  For `type = "HLIR"`, `"HLIRP"` or `"HLIRP (with host growth)"`, a data
  frame with one row per time step and the numeric columns `time`, `H`
  (healthy), `L` (latent), `I` (infectious) and `R` (removed), such as
  the output of
  [`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html) converted
  with [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).
  Other columns (for example `P`, the primary inoculum) are ignored.

  For `type = "MULTI"`, a list with the elements:

  - `n.yr`: number of seasons;

  - `results`: a list of `n.yr` data frames, one per season, each with
    the columns `time` (days since the start of the season), `H`, `L`,
    `I` and `R`;

  - `inoc`: numeric vector of length `n.yr` with the inoculum at the
    start of each season;

  - `n.max`: maximum amount of host tissue, used to set the upper limit
    of the y axis (`1.15 * n.max`).

- type:

  Model type, one of `"HLIR"`, `"HLIRP"`, `"HLIRP (with host growth)"`
  or `"MULTI"`.

- print:

  Logical. If `TRUE` (default) the plot(s) are printed, as in the book.
  Set to `FALSE` to only return them, for example to modify them with
  further ggplot2 layers.

## Value

Invisibly, a `ggplot` object for the HLIR-type models, or a list with
two `ggplot` objects, `progress` (diseased and total tissue) and
`inoculum` (log inoculum per season), for `type = "MULTI"`.

## Details

For the HLIR-type models, diseased tissue is computed as `D = L + I + R`
and drawn as a dotted line together with the four compartments. With
`type = "HLIRP (with host growth)"` the total amount of host tissue,
`T = H + L + I + R`, is added as a dotted blue line and the y axis is
labelled as an amount of host tissue instead of a number of hosts.

For `type = "MULTI"`, the seasons are placed one after the other on a
time axis in years (`((season - 1) * 365 + time) / 365`), with dashed
vertical lines at the end of each season. A second plot shows
`log(inoc)` at the start of each season.

## See also

[`theme_r4pde()`](https://emdelponte.github.io/r4pde/reference/theme_r4pde.md)

Other Disease modeling:
[`get_brdwgd()`](https://emdelponte.github.io/r4pde/reference/get_brdwgd.md),
[`get_era5()`](https://emdelponte.github.io/r4pde/reference/get_era5.md),
[`get_nasapower()`](https://emdelponte.github.io/r4pde/reference/get_nasapower.md),
[`windowpane()`](https://emdelponte.github.io/r4pde/reference/windowpane.md)

## Examples

``` r
# A simple HLIR model integrated with Euler steps
# (the book uses deSolve::ode())
hlir <- function(beta = 0.3, p = 5, i = 10, H0 = 99, I0 = 1,
                 days = 100, dt = 0.1) {
  n <- days / dt + 1
  out <- data.frame(time = seq(0, days, by = dt), H = H0, L = 0, I = I0, R = 0)
  for (k in 2:n) {
    s <- out[k - 1, ]
    inf <- beta * s$H * s$I / (s$H + s$L + s$I + s$R)
    out$H[k] <- s$H - dt * inf
    out$L[k] <- s$L + dt * (inf - s$L / p)
    out$I[k] <- s$I + dt * (s$L / p - s$I / i)
    out$R[k] <- s$R + dt * s$I / i
  }
  out
}
res <- hlir()
plot_model_results(res, type = "HLIR")


# Return the plot without printing, then modify it
p <- plot_model_results(res, type = "HLIR", print = FALSE)
p + ggplot2::labs(title = "HLIR model, beta = 0.3")


# Multi-seasonal epidemic: three seasons with increasing inoculum
multi <- list(
  n.yr = 3,
  results = lapply(1:3, function(k) hlir(I0 = k, days = 120, dt = 0.5)),
  inoc = c(1, 2, 3),
  n.max = 100
)
plot_model_results(multi, type = "MULTI")


```
