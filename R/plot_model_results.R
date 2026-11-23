#' Plot the output of temporal compartmental epidemic models
#'
#' @description
#' Helper used in the chapters on temporal compartmental models of the book
#' *R for Plant Disease Epidemiology* (R4PDE). It draws the dynamics of the
#' host compartments simulated by HLIR-type models, and the disease progress
#' and inoculum level across seasons of a multi-seasonal epidemic.
#'
#' @param data For `type = "HLIR"`, `"HLIRP"` or `"HLIRP (with host growth)"`,
#'   a data frame with one row per time step and the numeric columns `time`,
#'   `H` (healthy), `L` (latent), `I` (infectious) and `R` (removed), such as
#'   the output of `deSolve::ode()` converted with `as.data.frame()`. Other
#'   columns (for example `P`, the primary inoculum) are ignored.
#'
#'   For `type = "MULTI"`, a list with the elements:
#'   * `n.yr`: number of seasons;
#'   * `results`: a list of `n.yr` data frames, one per season, each with the
#'     columns `time` (days since the start of the season), `H`, `L`, `I` and `R`;
#'   * `inoc`: numeric vector of length `n.yr` with the inoculum at the start
#'     of each season;
#'   * `n.max`: maximum amount of host tissue, used to set the upper limit of
#'     the y axis (`1.15 * n.max`).
#' @param type Model type, one of `"HLIR"`, `"HLIRP"`,
#'   `"HLIRP (with host growth)"` or `"MULTI"`.
#' @param print Logical. If `TRUE` (default) the plot(s) are printed, as in the
#'   book. Set to `FALSE` to only return them, for example to modify them with
#'   further ggplot2 layers.
#'
#' @details
#' For the HLIR-type models, diseased tissue is computed as `D = L + I + R`
#' and drawn as a dotted line together with the four compartments. With
#' `type = "HLIRP (with host growth)"` the total amount of host tissue,
#' `T = H + L + I + R`, is added as a dotted blue line and the y axis is
#' labelled as an amount of host tissue instead of a number of hosts.
#'
#' For `type = "MULTI"`, the seasons are placed one after the other on a time
#' axis in years (`((season - 1) * 365 + time) / 365`), with dashed vertical
#' lines at the end of each season. A second plot shows `log(inoc)` at the
#' start of each season.
#'
#' @return Invisibly, a `ggplot` object for the HLIR-type models, or a list
#'   with two `ggplot` objects, `progress` (diseased and total tissue) and
#'   `inoculum` (log inoculum per season), for `type = "MULTI"`.
#'
#' @seealso [theme_r4pde()]
#'
#' @examples
#' # A simple HLIR model integrated with Euler steps
#' # (the book uses deSolve::ode())
#' hlir <- function(beta = 0.3, p = 5, i = 10, H0 = 99, I0 = 1,
#'                  days = 100, dt = 0.1) {
#'   n <- days / dt + 1
#'   out <- data.frame(time = seq(0, days, by = dt), H = H0, L = 0, I = I0, R = 0)
#'   for (k in 2:n) {
#'     s <- out[k - 1, ]
#'     inf <- beta * s$H * s$I / (s$H + s$L + s$I + s$R)
#'     out$H[k] <- s$H - dt * inf
#'     out$L[k] <- s$L + dt * (inf - s$L / p)
#'     out$I[k] <- s$I + dt * (s$L / p - s$I / i)
#'     out$R[k] <- s$R + dt * s$I / i
#'   }
#'   out
#' }
#' res <- hlir()
#' plot_model_results(res, type = "HLIR")
#'
#' # Return the plot without printing, then modify it
#' p <- plot_model_results(res, type = "HLIR", print = FALSE)
#' p + ggplot2::labs(title = "HLIR model, beta = 0.3")
#'
#' # Multi-seasonal epidemic: three seasons with increasing inoculum
#' multi <- list(
#'   n.yr = 3,
#'   results = lapply(1:3, function(k) hlir(I0 = k, days = 120, dt = 0.5)),
#'   inoc = c(1, 2, 3),
#'   n.max = 100
#' )
#' plot_model_results(multi, type = "MULTI")
#'
#' @export
#' @family Disease modeling
plot_model_results <- function(data,
                               type = c("HLIR", "HLIRP",
                                        "HLIRP (with host growth)", "MULTI"),
                               print = TRUE) {
  type <- match.arg(type)

  if (type == "MULTI") {
    out <- .plot_multi_season(data)
    if (isTRUE(print)) {
      print(out$progress)
      print(out$inoculum)
    }
    return(invisible(out))
  }

  p <- .plot_hlir(data, type)
  if (isTRUE(print)) print(p)
  invisible(p)
}

# Check that `x` is a data frame with the compartment columns
.check_compartments <- function(x, what = "`data`") {
  if (!is.data.frame(x)) {
    stop(what, " must be a data frame (use as.data.frame() on the output of deSolve::ode()).",
         call. = FALSE)
  }
  need <- c("time", "H", "L", "I", "R")
  miss <- setdiff(need, names(x))
  if (length(miss) > 0) {
    stop(what, " must have the columns ", paste(need, collapse = ", "),
         "; missing: ", paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

.plot_hlir <- function(data, type) {
  .check_compartments(data)
  plot_df <- as.data.frame(data)
  plot_df$D <- plot_df$L + plot_df$I + plot_df$R

  cols <- c(H = "darkolivegreen", L = "salmon", I = "firebrick",
            R = "grey40", D = "black")
  ltys <- c(H = "solid", L = "solid", I = "solid", R = "solid", D = "dotted")
  labels <- c(H = "(H)ealthy", L = "(L)atent", I = "(I)nfectious",
              R = "(R)emoved", D = "(D)iseased")
  y_lab <- "Number of hosts"

  if (type == "HLIRP (with host growth)") {
    plot_df$T_tissue <- plot_df$H + plot_df$L + plot_df$I + plot_df$R
    cols <- c(cols, T_tissue = "dodgerblue")
    ltys <- c(ltys, T_tissue = "dotted")
    labels <- c(labels, T_tissue = "(T)otal")
    y_lab <- "Amount of host tissue"
  }
  keep <- names(cols)

  plot_long <- plot_df |>
    dplyr::select(dplyr::all_of(c("time", keep))) |>
    tidyr::pivot_longer(cols = -"time", names_to = "Compartment",
                        values_to = "Value") |>
    dplyr::mutate(Compartment = factor(.data$Compartment, levels = keep))

  ggplot2::ggplot(plot_long,
                  ggplot2::aes(x = .data$time, y = .data$Value,
                               color = .data$Compartment,
                               linetype = .data$Compartment)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = cols, labels = labels) +
    ggplot2::scale_linetype_manual(values = ltys, labels = labels) +
    ggplot2::labs(x = "Time (days)", y = y_lab,
                  title = paste(type, "model"),
                  color = NULL, linetype = NULL) +
    theme_r4pde() +
    .legend_inside(c(0.78, 0.5)) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 8))
}

# Legend inside the panel; ggplot2 >= 3.5.0 deprecates numeric legend.position
.legend_inside <- function(pos) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    ggplot2::theme(legend.position = "inside", legend.position.inside = pos)
  } else {
    ggplot2::theme(legend.position = pos)
  }
}

.plot_multi_season <- function(data) {
  need <- c("n.yr", "results", "inoc", "n.max")
  if (!is.list(data) || is.data.frame(data) || !all(need %in% names(data))) {
    stop("For type = \"MULTI\", `data` must be a list with the elements ",
         paste(need, collapse = ", "), ".", call. = FALSE)
  }
  n_yr <- data$n.yr
  if (length(data$results) < n_yr) {
    stop("`data$results` has ", length(data$results), " season(s) but `data$n.yr` is ",
         n_yr, ".", call. = FALSE)
  }
  if (length(data$inoc) < n_yr) {
    stop("`data$inoc` must have one value per season (", n_yr, ").", call. = FALSE)
  }

  breaks <- numeric(n_yr)
  multi_list <- vector("list", n_yr)
  for (i in seq_len(n_yr)) {
    res <- data$results[[i]]
    .check_compartments(res, paste0("`data$results[[", i, "]]`"))
    t_years <- ((i - 1) * 365 + res$time) / 365
    multi_list[[i]] <- data.frame(
      time = t_years,
      D = res$L + res$I + res$R,
      T_tissue = res$H + res$L + res$I + res$R,
      Season = factor(i)
    )
    breaks[i] <- max(t_years)
  }

  multi_long <- dplyr::bind_rows(multi_list) |>
    tidyr::pivot_longer(cols = c("D", "T_tissue"), names_to = "Variable",
                        values_to = "Value")

  p1 <- ggplot2::ggplot(multi_long,
                        ggplot2::aes(x = .data$time, y = .data$Value,
                                     color = .data$Variable,
                                     group = interaction(.data$Season, .data$Variable))) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = breaks, color = "grey",
                        linetype = "dashed", linewidth = 0.8) +
    ggplot2::scale_color_manual(
      values = c(D = "black", T_tissue = "dodgerblue"),
      labels = c(D = "(D)iseased", T_tissue = "(T)otal")
    ) +
    ggplot2::labs(x = "Time (years)", y = "Amount of host tissue",
                  title = "Multi-seasonal epidemic. Disease progress",
                  color = NULL) +
    ggplot2::ylim(0, 1.15 * data$n.max) +
    theme_r4pde() +
    ggplot2::theme(legend.position = "top")

  inoc_df <- data.frame(Year = seq_len(n_yr), log_Inoc = log(data$inoc[seq_len(n_yr)]))
  p2 <- ggplot2::ggplot(inoc_df, ggplot2::aes(x = .data$Year, y = .data$log_Inoc)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = "Time (years)",
                  y = "log(Inoculum) at start of season",
                  title = "Multi-seasonal epidemic. log(Inoculum) level") +
    ggplot2::scale_x_continuous(breaks = seq_len(n_yr)) +
    theme_r4pde()

  list(progress = p1, inoculum = p2)
}
