# Load necessary libraries
library(cowplot)
library(ggplot2)

#' Custom ggplot2 theme based on cowplot::theme_half_open
#'
#' This function creates a new ggplot2 theme by modifying the cowplot::theme_half_open theme.
#' It sets a custom font size and changes the panel background color to gray96.
#'
#' @param font_size The base font size. Default is 16.
#'
#' @return A ggplot2 theme object.
#' @export

theme_r4pde <- function(font_size = 16) {
  base_theme <- cowplot::theme_half_open(font_size = font_size)
  modified_theme <- base_theme +
    theme(panel.background = element_rect(fill = "gray96"),
          strip.background = element_rect(colour = "white", fill = "white"))
  return(modified_theme)
}
