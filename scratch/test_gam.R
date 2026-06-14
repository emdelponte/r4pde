library(dplyr)
library(ggplot2)
library(mgcv)
library(tidyr)

source("R/functional_curves.R")
source("R/functional_dsp.R")

dat <- simulate_dsp_data()
dat_gam <- dat |> rename(y = severity, trt = treatment) |> mutate(trt = as.factor(trt))

# Try original gam
res <- functional_curves(dat_gam, time="time", response="y", treatment="trt", response_scale="percent")
ggsave("test_original.png", res$plot_mean)

# Let's override the function to test without global smooth
functional_curves_modified <- functional_curves
body(functional_curves_modified) <- bquote({
  .(body(functional_curves_modified))
})

# Actually, I'll just change the text
code <- readLines("R/functional_curves.R")
code <- gsub('sprintf("s(%s, k=%s)", .time, k_smooth_eff),', '', code, fixed=TRUE)
code <- gsub('sprintf("s(%s, by=%s, k=%s)", .time, .trt, k_trt_eff)', 'sprintf("s(%s, by=%s, k=%s)", .time, .trt, k_smooth_eff)', code, fixed=TRUE)
writeLines(code, "scratch/functional_curves_mod.R")

source("scratch/functional_curves_mod.R")
res_mod <- functional_curves(dat_gam, time="time", response="y", treatment="trt", response_scale="percent")
ggsave("test_mod.png", res_mod$plot_mean)
