#' Normalized functional instability from compare_curves output
#'
#' Computes normalized functional instability (NFI) for each treatment
#' based on genotype-by-environment predicted curves extracted from a
#' `compare_curves()` object. Optionally, instability can be decomposed
#' into spatial and temporal components if the environment identifier can
#' be split into location and year.
#'
#' The overall instability metric is defined as the mean integrated squared
#' deviation of each genotype-by-environment curve from the genotype-specific
#' mean curve across environments, normalized by the integrated squared mean
#' curve.
#'
#' @param x An object returned by [compare_curves()].
#' @param n_time Number of points in the prediction grid over the time domain.
#' @param env_sep Optional separator used to split `env` into spatial and
#'   temporal components, for example `"_"` or `"-"`. If `NULL`, only the
#'   overall NFI is returned.
#' @param env_names A character vector of length 2 giving the names of the
#'   components obtained after splitting `env`. Default is
#'   `c("location", "year")`.
#' @param return_curves Logical. If `TRUE`, the predicted genotype-by-environment
#'   curves are also returned.
#'
#' @return If `return_curves = FALSE`, a tibble with one row per genotype and
#'   columns:
#'   \describe{
#'     \item{geno}{Treatment or genotype identifier.}
#'     \item{n_env}{Number of environments used in the calculation.}
#'     \item{FI}{Absolute functional instability.}
#'     \item{mean_energy}{Integrated squared mean epidemic curve.}
#'     \item{nFI}{Normalized functional instability. Lower values indicate
#'       greater stability.}
#'   }
#'
#'   If `env_sep` is provided, the returned tibble also includes:
#'   \describe{
#'     \item{FI_space, mean_energy_space, nFI_space}{Spatial component of
#'       instability.}
#'     \item{FI_time, mean_energy_time, nFI_time}{Temporal component of
#'       instability.}
#'   }
#'
#'   If `return_curves = TRUE`, a list is returned with two elements:
#'   \describe{
#'     \item{metrics}{The tibble described above.}
#'     \item{curves}{A tibble of predicted genotype-by-environment curves with
#'       columns `geno`, `env`, `time`, `mu`, and `eta`.}
#'   }
#'
#' @details
#' Let \eqn{f_{ge}(t)} denote the predicted epidemic trajectory of genotype
#' \eqn{g} in environment \eqn{e}, and let \eqn{\bar f_g(t)} denote the mean
#' trajectory of genotype \eqn{g} across environments. Functional instability
#' is computed as:
#'
#' \deqn{
#' FI_g = \frac{1}{E_g}\sum_{e=1}^{E_g} \int_T \left(f_{ge}(t)-\bar f_g(t)\right)^2 dt
#' }
#'
#' and normalized as:
#'
#' \deqn{
#' nFI_g = \frac{FI_g}{\int_T \bar f_g(t)^2 dt}
#' }
#'
#' Numerical integration is performed with the trapezoidal rule on a regular
#' prediction grid over the observed time domain.
#'
#' @seealso [compare_curves()]
#'
#' @examples
#' \dontrun{
#' m1 <- r4pde::compare_curves(
#'   data = dat_ready,
#'   time = "time",
#'   response = "y",
#'   treatment = "geno",
#'   environment = "env",
#'   cluster_k = 4
#' )
#'
#' # Overall instability
#' res <- functional_instability(m1)
#' res
#'
#' # Overall + spatial and temporal components
#' # Assuming env has values such as "PF_2021"
#' res2 <- functional_instability(m1, env_sep = "_")
#' res2
#'
#' # Return curves used in the calculations
#' out <- functional_instability(m1, env_sep = "_", return_curves = TRUE)
#' out$metrics
#' head(out$curves)
#' }
#'
#' @export
functional_instability <- function(x,
                n_time = 200,
                env_sep = NULL,
                env_names = c("location", "year"),
                return_curves = FALSE) {

  if (is.null(names(x)) || !all(c("gam", "data") %in% names(x))) {
    stop("`x` must be a valid `compare_curves()` object containing `gam` and `data`.",
         call. = FALSE)
  }

  dat <- x$data

  trt_col <- if (!is.null(x$vars$treatment)) x$vars$treatment else "geno"
  env_col <- if (!is.null(x$vars$environment)) x$vars$environment else "env"
  time_col <- if (!is.null(x$vars$time)) x$vars$time else "time"

  if (!all(c(trt_col, env_col, time_col) %in% names(dat))) {
    stop(sprintf("`x$data` must contain columns `%s`, `%s`, and `%s`.", trt_col, env_col, time_col), call. = FALSE)
  }

 trapz_vec <- function(x, y) {
  ok <- !(is.na(x) | is.na(y))
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 2) return(NA_real_)

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

  predict_geno_env_curves <- function(dat, mod, trt_c, env_c, time_c, n_time = 200) {
    pairs <- dat %>%
      dplyr::distinct(.data[[trt_c]], .data[[env_c]])

    time_seq <- seq(
      min(dat[[time_c]], na.rm = TRUE),
      max(dat[[time_c]], na.rm = TRUE),
      length.out = n_time
    )

    time_df <- tibble::tibble(time_seq)
    names(time_df) <- time_c

    newdata <- tidyr::expand_grid(pairs, time_df)

    mu  <- stats::predict(mod, newdata = newdata, type = "response")
    eta <- stats::predict(mod, newdata = newdata, type = "link")

    newdata <- newdata %>%
      dplyr::mutate(mu = as.numeric(mu), eta = as.numeric(eta))

    names(newdata)[names(newdata) == trt_c] <- "geno"
    names(newdata)[names(newdata) == env_c] <- "env"
    names(newdata)[names(newdata) == time_c] <- "time"

    newdata
  }

  calc_nfi_general <- function(curves) {
    mean_curve <- curves %>%
      dplyr::group_by(geno, time) %>%
      dplyr::summarise(mean_mu = mean(mu, na.rm = TRUE), .groups = "drop")

    fi_env <- curves %>%
      dplyr::left_join(mean_curve, by = c("geno", "time")) %>%
      dplyr::mutate(sq_dev = (mu - mean_mu)^2) %>%
      dplyr::group_by(geno, env) %>%
      dplyr::summarise(
        int_sq_dev = trapz_vec(time, sq_dev),
        .groups = "drop"
      )

    fi_geno <- fi_env %>%
      dplyr::group_by(geno) %>%
      dplyr::summarise(
        n_env = dplyr::n(),
        FI = mean(int_sq_dev, na.rm = TRUE),
        .groups = "drop"
      )

    mean_energy <- mean_curve %>%
      dplyr::group_by(geno) %>%
      dplyr::summarise(
        mean_energy = trapz_vec(time, mean_mu^2),
        .groups = "drop"
      )

    fi_geno %>%
      dplyr::left_join(mean_energy, by = "geno") %>%
      dplyr::mutate(nFI = FI / mean_energy) %>%
      dplyr::arrange(nFI)
  }

  calc_nfi_space <- function(curves) {
    mean_curve <- curves %>%
      dplyr::group_by(geno, year, time) %>%
      dplyr::summarise(mean_mu = mean(mu, na.rm = TRUE), .groups = "drop")

    fi_units <- curves %>%
      dplyr::left_join(mean_curve, by = c("geno", "year", "time")) %>%
      dplyr::mutate(sq_dev = (mu - mean_mu)^2) %>%
      dplyr::group_by(geno, year, location) %>%
      dplyr::summarise(
        int_sq_dev = trapz_vec(time, sq_dev),
        .groups = "drop"
      )

    fi_geno <- fi_units %>%
      dplyr::group_by(geno) %>%
      dplyr::summarise(
        FI_space = mean(int_sq_dev, na.rm = TRUE),
        .groups = "drop"
      )

    mean_energy <- mean_curve %>%
      dplyr::group_by(geno, year) %>%
      dplyr::summarise(
        mean_energy_space = trapz_vec(time, mean_mu^2),
        .groups = "drop"
      ) %>%
      dplyr::group_by(geno) %>%
      dplyr::summarise(
        mean_energy_space = mean(mean_energy_space, na.rm = TRUE),
        .groups = "drop"
      )

    fi_geno %>%
      dplyr::left_join(mean_energy, by = "geno") %>%
      dplyr::mutate(nFI_space = FI_space / mean_energy_space)
  }

  calc_nfi_time <- function(curves) {
    mean_curve <- curves %>%
      dplyr::group_by(geno, location, time) %>%
      dplyr::summarise(mean_mu = mean(mu, na.rm = TRUE), .groups = "drop")

    fi_units <- curves %>%
      dplyr::left_join(mean_curve, by = c("geno", "location", "time")) %>%
      dplyr::mutate(sq_dev = (mu - mean_mu)^2) %>%
      dplyr::group_by(geno, location, year) %>%
      dplyr::summarise(
        int_sq_dev = trapz_vec(time, sq_dev),
        .groups = "drop"
      )

    fi_geno <- fi_units %>%
      dplyr::group_by(geno) %>%
      dplyr::summarise(
        FI_time = mean(int_sq_dev, na.rm = TRUE),
        .groups = "drop"
      )

    mean_energy <- mean_curve %>%
      dplyr::group_by(geno, location) %>%
      dplyr::summarise(
        mean_energy_time = trapz_vec(time, mean_mu^2),
        .groups = "drop"
      ) %>%
      dplyr::group_by(geno) %>%
      dplyr::summarise(
        mean_energy_time = mean(mean_energy_time, na.rm = TRUE),
        .groups = "drop"
      )

    fi_geno %>%
      dplyr::left_join(mean_energy, by = "geno") %>%
      dplyr::mutate(nFI_time = FI_time / pmax(mean_energy_time, 1e-8)) %>%
      dplyr::mutate(nFI_time = pmax(nFI_time, 0))
  }

  curves <- predict_geno_env_curves(dat, x$gam, trt_col, env_col, time_col, n_time = n_time)
  out <- calc_nfi_general(curves)

  if (!is.null(env_sep)) {
    if (length(env_names) != 2) {
      stop("`env_names` must have length 2.", call. = FALSE)
    }

    curves2 <- curves %>%
      tidyr::separate(
        env,
        into = env_names,
        sep = env_sep,
        remove = FALSE
      )

    names(curves2)[names(curves2) == env_names[1]] <- "location"
    names(curves2)[names(curves2) == env_names[2]] <- "year"

    out <- out %>%
      dplyr::left_join(calc_nfi_space(curves2), by = "geno") %>%
      dplyr::left_join(calc_nfi_time(curves2), by = "geno")
  }

  class(out) <- c("fi_tbl", class(out))

  if ("geno" %in% names(out)) {
    names(out)[names(out) == "geno"] <- trt_col
  }

  if (return_curves) {
    names(curves)[names(curves) == "geno"] <- trt_col
    names(curves)[names(curves) == "env"] <- env_col
    names(curves)[names(curves) == "time"] <- time_col
    return(list(metrics = out, curves = curves))
  }

  out
}
