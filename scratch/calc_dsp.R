t <- seq(0, 60, by = 1)

control_curve <- function(t) {
  78 / (1 + exp(-0.13 * (t - 35)))
}
y_control <- control_curve(t)

y_early <- 76 / (1 + exp(-0.2 * (t - 48)))
y_late <- 25 / (1 + exp(-0.15 * (t - 25)))
y_persist <- 35 / (1 + exp(-0.08 * (t - 35)))

get_area <- function(y) sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2)

cat("AUDPC Early:", get_area(y_early), "\n")
cat("AUDPC Late:", get_area(y_late), "\n")
cat("AUDPC Persist:", get_area(y_persist), "\n")

dsp_e <- y_control - y_early
dsp_l <- y_control - y_late
dsp_p <- y_control - y_persist

cat("Max DSP Early at t=", t[which.max(dsp_e)], " (", max(dsp_e), ")\n")
cat("Max DSP Late at t=", t[which.max(dsp_l)], " (", max(dsp_l), ")\n")
cat("Max DSP Persist at t=", t[which.max(dsp_p)], " (", max(dsp_p), ")\n")
