# Instantaneous Disease Progress Rates from Functional Curves

Estimates instantaneous rates of plant disease progress (\\S'(t) =
dS(t)/dt\\) and associated uncertainty from epidemic trajectories fitted
by
[`functional_curves`](https://emdelponte.github.io/r4pde/reference/functional_curves.md).
Derivatives are computed using the GAM linear predictor matrix
(`lpmatrix`) with finite differences, supporting both the original
response scale (e.g., severity proportion/percentage points per day) and
the link (linear predictor) scale. Also computes epidemic rate
phenotypes including maximum rate (\\r\_{max}\\), time of maximum rate
(\\t\_{r\_{max}}\\), growth duration, cumulative positive growth, and
distinct growth windows.

## Usage

``` r
functional_rate(
  object,
  scale = c("response", "link"),
  method = c("central"),
  n_grid = 200,
  eps = NULL,
  interval = 0.95,
  threshold = 0,
  positive_only = TRUE,
  uncertainty = TRUE,
  growth_criterion = c("rate", "significant"),
  newdata = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `"functional_curves"`.

- scale:

  Character specifying the derivative scale: `"response"` (default;
  e.g., severity change per unit time) or `"link"` (linear predictor
  scale).

- method:

  Character specifying the finite difference method. Currently
  `"central"` (central difference in the interior with one-sided
  differences at temporal boundaries).

- n_grid:

  Integer; number of points in the regular temporal evaluation grid per
  group. Default is 200.

- eps:

  Numeric step size used in finite differences. If `NULL` (default),
  automatically chosen as \\10^{-4} \times (t\_{max} - t\_{min})\\.

- interval:

  Numeric confidence level for intervals. Default is 0.95.

- threshold:

  Numeric minimum rate threshold considered biologically relevant.
  Default is 0.

- positive_only:

  Logical; whether to restrict summarized growth descriptors to positive
  rates. Default is `TRUE`. Note that negative rates are always
  preserved in the detailed rate table (`data`).

- uncertainty:

  Logical; whether to compute standard errors and confidence intervals
  using the GAM covariance matrix. Default is `TRUE`.

- growth_criterion:

  Character specifying how active epidemic growth is defined: `"rate"`
  (default; `rate > threshold`) or `"significant"`
  (`rate_lower > threshold`, requiring `uncertainty = TRUE`).

- newdata:

  Optional user-supplied data frame with custom evaluation points. Must
  contain the time and treatment variables.

- ...:

  Additional arguments passed to methods.

## Value

An object of class `"functional_rate"` containing:

- `data`:

  A tibble of full evaluation grid with columns including treatment,
  time, `fitted`, `fitted_se`, `fitted_lower`, `fitted_upper`, `rate`,
  `rate_se`, `rate_lower`, `rate_upper`, `significant_increase`,
  `above_threshold`, and `significant_above_threshold`.

- `summary`:

  A tibble of curve-level rate phenotypes: `r_max`, `t_r_max`, `r_mean`,
  `r_mean_positive`, `t_growth_start`, `t_growth_end`,
  `growth_duration`, `cumulative_positive_growth`, and
  `n_growth_windows`.

- `model`:

  The underlying GAM model object.

- `call`:

  The matched call.

- `settings`:

  List of analysis settings used.

- `vars`:

  List of variable names identified in the model.

- `family`:

  Character string naming the GAM family.

## Details

**Mathematical Formulation:** Let \\S(t)\\ denote the fitted disease
trajectory and \\\eta(t) = X(t)\beta\\ the linear predictor from the GAM
model, with \\\mu(t) = g^{-1}(\eta(t))\\.

The derivative on the link scale is estimated by: \$\$D = \frac{X(t +
\epsilon) - X(t - \epsilon)}{2\epsilon}\$\$ \$\$\eta'(t) = D \beta\$\$
with variance \\\operatorname{Var}(\eta'(t)) = \operatorname{diag}(D V_p
D^T)\\, where \\V_p\\ is the Bayesian posterior covariance matrix of the
GAM parameters.

On the response scale (`scale = "response"`), the derivative of the mean
trajectory is: \$\$S'(t) = \frac{\mu(\eta(t + \epsilon)) - \mu(\eta(t -
\epsilon))}{2\epsilon}\$\$ The approximate gradient vector \\G\\ with
respect to \\\beta\\ is: \$\$G = \frac{\mu'(\eta(t + \epsilon)) X(t +
\epsilon) - \mu'(\eta(t - \epsilon)) X(t - \epsilon)}{2\epsilon}\$\$
where \\\mu'(\eta) = d\mu/d\eta\\ is given by `family$mu.eta()`. The
uncertainty is then: \$\$\operatorname{Var}(S'(t)) =
\operatorname{diag}(G V_p G^T)\$\$

**Boundary Handling:** At the observed boundaries \\t\_{min}\\ and
\\t\_{max}\\, one-sided forward and backward differences are used to
prevent unintended extrapolation beyond the observed domain.

**Time-Invariant Terms and Random Effects:** Intercepts, main treatment
effects, and terms constant with respect to time cancel analytically in
\\X(t + \epsilon) - X(t - \epsilon)\\. Environmental and experimental
unit random effects are excluded by default to yield population-adjusted
treatment rates, consistent with
[`functional_curves`](https://emdelponte.github.io/r4pde/reference/functional_curves.md).

**No-Epidemic and Flat Curves:** Curves with all observed response
values equal to zero (or estimated rates below numerical tolerance
\\10^{-6}\\) explicitly return \\r\_{max} = 0\\, \\t\_{r\_{max}} =
\text{NA}\\, \\\text{growth\\duration} = 0\\, and
\\\text{cumulative\\positive\\growth} = 0\\.

## References

Madden, L. V., Hughes, G., & van den Bosch, F. (2007). The Study of
Plant Disease Epidemics. APS Press, St. Paul, MN.

Wood, S. N. (2017). Generalized Additive Models: An Introduction with R
(2nd ed.). Chapman and Hall/CRC.

Simpson, G. L. (2018). Modelling Palaeoecological Time Series Using
Generalized Additive Models. Frontiers in Ecology and Evolution, 6, 149.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
df <- data.frame(
  time = rep(1:30, 3),
  treatment = rep(c("Early", "Late", "None"), each = 30),
  y = c(
    1 / (1 + exp(-(1:30 - 10) * 0.4)),  # Early
    1 / (1 + exp(-(1:30 - 20) * 0.4)),  # Late
    rep(0, 30)                          # None
  )
)

curves <- functional_curves(
  data = df, time = "time", response = "y", treatment = "treatment",
  min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
)

rates <- functional_rate(curves)
print(rates)
summary(rates)
plot(rates)
} # }
```
