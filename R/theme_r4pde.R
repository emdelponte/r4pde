#' Custom ggplot2 theme based on cowplot::theme_half_open
#'
#' This function creates a new ggplot2 theme by modifying the cowplot::theme_half_open theme.
#' It sets a custom font size and changes the panel background color to gray96.
#' If cowplot is not installed, a close equivalent built on
#' \code{ggplot2::theme_classic()} is used.
#'
#' @param font_size The base font size. Default is 16.
#'
#' @return A ggplot2 theme object.
#' @export
#' @family Miscellaneous

theme_r4pde <- function(font_size = 16) {
  base_theme <- if (requireNamespace("cowplot", quietly = TRUE)) {
    cowplot::theme_half_open(font_size = font_size)
  } else {
    ggplot2::theme_classic(base_size = font_size) +
      ggplot2::theme(
        axis.text = ggplot2::element_text(colour = "black", size = ggplot2::rel(0.86)),
        axis.ticks = ggplot2::element_line(colour = "black"),
        axis.line = ggplot2::element_line(colour = "black")
      )
  }
  modified_theme <- base_theme +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "gray96"),
                   strip.background = ggplot2::element_rect(colour = "white", fill = "white"))
  return(modified_theme)
}
