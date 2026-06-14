time_seq <- seq(0, 60, by = 5)
control_curve <- function(t) {
  78 / (1 + exp(-0.13 * (t - 35)))
}
y_control <- control_curve(time_seq)

trapz <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

audpc_control <- trapz(time_seq, y_control)
target_audpc <- 0.30 * audpc_control

cat("Control AUDPC:", audpc_control, "\n")
cat("Target AUDPC (30%):", target_audpc, "\n")

# Objective function to find the right amplitude
obj_early <- function(amp) {
  dsp <- amp * exp(-((time_seq - 30)^2) / (2 * 12^2))
  sev <- pmax(0, y_control - dsp)
  sev <- cummax(sev)
  abs(trapz(time_seq, sev) - target_audpc)
}

obj_late <- function(amp) {
  dsp <- amp * exp(-((time_seq - 58)^2) / (2 * 15^2))
  sev <- pmax(0, y_control - dsp)
  sev <- cummax(sev)
  abs(trapz(time_seq, sev) - target_audpc)
}

obj_persist <- function(amp) {
  dsp <- amp * plogis((time_seq - 20) / 3) * (1 - plogis((time_seq - 65) / 5))
  sev <- pmax(0, y_control - dsp)
  sev <- cummax(sev)
  abs(trapz(time_seq, sev) - target_audpc)
}

opt_early <- optimize(obj_early, interval = c(0, 500))
opt_late <- optimize(obj_late, interval = c(0, 500))
opt_persist <- optimize(obj_persist, interval = c(0, 500))

cat("Optimal amp Early:", opt_early$minimum, "\n")
cat("Optimal amp Late:", opt_late$minimum, "\n")
cat("Optimal amp Persist:", opt_persist$minimum, "\n")
