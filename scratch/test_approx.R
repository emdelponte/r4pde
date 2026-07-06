library(dplyr)
dat <- tibble(trt = c("A", "A", "A", "A"), time = 1:4, y = c(2,3,4,5))
dat |> 
  group_by(trt) |>
  reframe(
    y = stats::approx(x = .data$time, y = .data$y, xout = c(1,2,3))$y,
    time = c(1,2,3)
  )
