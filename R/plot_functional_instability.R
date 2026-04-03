#' Plot method for functional instability results
#'
#' @param object Output from [functional_instability()].
#' @param type One of `"overall"`, `"space"`, or `"time"`.
#' @param ... Further arguments passed to methods.
#'
#' @return A ggplot object.
#' @export
plot_functional_instability <- function(object, type = c("overall", "space", "time"), ...) {
  type <- match.arg(type)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_functional_instability().", call. = FALSE)
  }

  dat <- tibble::as_tibble(object)

  if (type == "overall") {
    p <- dat |>
      dplyr::mutate(geno = stats::reorder(geno, nFI)) |>
      ggplot2::ggplot(ggplot2::aes(geno, nFI)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = "Genotype",
        y = "Normalized functional instability (nFI)"
      ) +
      ggplot2::theme_minimal()
    return(p)
  }

  if (type == "space") {
    if (!"nFI_space" %in% names(dat)) {
      stop("`nFI_space` not found. Re-run `functional_instability()` with `env_sep`.", call. = FALSE)
    }

    p <- dat |>
      dplyr::mutate(geno = stats::reorder(geno, nFI_space)) |>
      ggplot2::ggplot(ggplot2::aes(geno, nFI_space)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = "Genotype",
        y = "Spatial normalized functional instability (nFI_space)"
      ) +
      ggplot2::theme_minimal()
    return(p)
  }

  if (!"nFI_time" %in% names(dat)) {
    stop("`nFI_time` not found. Re-run `functional_instability()` with `env_sep`.", call. = FALSE)
  }

  dat |>
    dplyr::mutate(geno = stats::reorder(geno, nFI_time)) |>
    ggplot2::ggplot(ggplot2::aes(geno, nFI_time)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "Genotype",
      y = "Temporal normalized functional instability (nFI_time)"
    ) +
    ggplot2::theme_minimal()
}
