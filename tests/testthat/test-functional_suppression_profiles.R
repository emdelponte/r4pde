test_that("functional_suppression_profiles works correctly", {
  sim_dat <- tibble::tibble(
    treatment = rep(c("Control", "A", "B", "C"), each = 6),
    time = rep(seq(0, 25, by = 5), times = 4),
    severity = c(
      c(5, 10, 20, 35, 50, 65),
      c(3, 5, 10, 18, 30, 40),
      c(4, 6, 12, 22, 35, 45),
      c(2, 4, 8, 14, 22, 32)
    )
  )

  fsp <- functional_suppression_profiles(
    data = sim_dat,
    reference = "Control",
    time = "time",
    response = "severity",
    treatment = "treatment",
    threshold = 5,
    k = 2
  )

  # Check outputs
  expect_s3_class(fsp, "functional_suppression_profiles")
  expect_true("classification" %in% names(fsp))
  expect_true(all(c("treatment", "profile", "mean_rank") %in% names(fsp$classification)))
  
  # Check mean_rank order
  expect_false(is.unsorted(fsp$classification$mean_rank))
  
  # Check profile names
  expect_true(all(grepl("^Profile ", as.character(fsp$classification$profile))))
  
  # 1. Check object class
  expect_s3_class(fsp, "functional_suppression_profiles")
  
  # 2. Check components
  expect_true(!is.null(fsp$contrast))
  expect_true(!is.null(fsp$dsp_curves))
  expect_true(!is.null(fsp$distances))
  expect_true(!is.null(fsp$hclust))
  expect_true(!is.null(fsp$clusters))
  expect_true(!is.null(fsp$summary))
  expect_true(!is.null(fsp$ranking))
  
  # Check independent ranking and summary
  expect_true("protected_area" %in% names(fsp$summary))
  expect_false("protected_area" %in% names(fsp$ranking))
  expect_true("rank_protected_area" %in% names(fsp$ranking))
  expect_true("mean_rank" %in% names(fsp$ranking))
  expect_false("average_rank" %in% names(fsp$ranking))
  
  # Check distances uses smoothed curves
  expect_identical(fsp$distances$functional_curves, fsp$dsp_curves)
  
  # 3. Check cluster table columns
  expect_true(all(c("treatment", "cluster") %in% names(fsp$clusters)))
  
  # 4. Check print method returns invisible object
  expect_invisible(print(fsp))
  
  # 5. Check summary method returns expected components
  summ <- summary(fsp)
  expect_s3_class(summ, "summary.functional_suppression_profiles")
  expect_true(all(c("classification", "profile_summary", "ranking", "metrics", "silhouette") %in% names(summ)))
  
  # 6. Check plot method returns a ggplot object
  p1 <- plot(fsp, type = "dendrogram")
  expect_s3_class(p1, "ggplot")
  
  p2 <- plot(fsp, type = "profiles")
  expect_s3_class(p2, "ggplot")
  expect_true("mu" %in% names(p2$data)) # smoothed curves
  
  p3 <- plot(fsp, type = "heatmap")
  expect_s3_class(p3, "ggplot")
  expect_true("rank" %in% names(p3$data)) # only ranks
  
  p4 <- plot(fsp, type = "rank")
  expect_s3_class(p4, "ggplot")
  expect_true("mean_rank" %in% names(p4$data))
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_all <- plot(fsp, type = "all")
    expect_s3_class(p_all, "patchwork")
  }
  
  # 7. missing metric handling works
  fsp_missing_metric <- functional_suppression_profiles(
    data = sim_dat,
    reference = "Control",
    time = "time",
    response = "severity",
    treatment = "treatment",
    metrics = c("protected_area", "non_existent_metric"),
    k = 2
  )
  expect_true("rank_protected_area" %in% names(fsp_missing_metric$ranking))
  expect_false("rank_non_existent_metric" %in% names(fsp_missing_metric$ranking))
  
  # 8. non-numeric treatment labels work
  sim_dat_char <- sim_dat
  sim_dat_char$treatment <- paste0("Trt_", sim_dat_char$treatment)
  fsp_char <- functional_suppression_profiles(
    data = sim_dat_char,
    reference = "Trt_Control",
    time = "time",
    response = "severity",
    treatment = "treatment",
    threshold = 5,
    k = 2
  )
  expect_true(all(grepl("Trt_", fsp_char$clusters$treatment)))
})

test_that("functional_suppression_profiles supports multiple environments", {
  # Create data with two environments
  sim_dat1 <- tibble::tibble(
    treatment = rep(c("Control", "A", "B", "C"), each = 6),
    time = rep(seq(0, 25, by = 5), times = 4),
    severity = c(
      c(5, 10, 20, 35, 50, 65),
      c(3, 5, 10, 18, 30, 40),
      c(4, 6, 12, 22, 35, 45),
      c(2, 4, 8, 14, 22, 32)
    ),
    env = "Env1"
  )
  
  sim_dat2 <- tibble::tibble(
    treatment = rep(c("Control", "A", "B", "C"), each = 6),
    time = rep(seq(0, 25, by = 5), times = 4),
    severity = c(
      c(6, 12, 24, 38, 54, 70),
      c(4, 6, 12, 20, 32, 44),
      c(5, 7, 14, 24, 38, 50),
      c(3, 5, 9, 16, 25, 36)
    ),
    env = "Env2"
  )
  
  sim_dat_multi <- dplyr::bind_rows(sim_dat1, sim_dat2)
  
  # Test Joint Mode
  fsp_joint <- functional_suppression_profiles(
    data = sim_dat_multi,
    reference = "Control",
    time = "time",
    response = "severity",
    treatment = "treatment",
    environment = "env",
    by_environment = FALSE,
    threshold = 5,
    k = 2
  )
  
  expect_s3_class(fsp_joint, "functional_suppression_profiles")
  expect_equal(nrow(fsp_joint$classification), 3) # A, B, C (reference is excluded)
  expect_true("mean_rank" %in% names(fsp_joint$classification))
  expect_false(any(duplicated(fsp_joint$classification$treatment)))
  
  # Test Separate Mode
  fsp_sep <- functional_suppression_profiles(
    data = sim_dat_multi,
    reference = "Control",
    time = "time",
    response = "severity",
    treatment = "treatment",
    environment = "env",
    by_environment = TRUE,
    threshold = 5,
    k = 2
  )
  
  expect_s3_class(fsp_sep, "functional_suppression_profiles_list")
  expect_length(fsp_sep, 2)
  expect_identical(names(fsp_sep), c("Env1", "Env2"))
  expect_s3_class(fsp_sep$Env1, "functional_suppression_profiles")
  expect_s3_class(fsp_sep$Env2, "functional_suppression_profiles")
  
  # Test print method for the list
  expect_output(print(fsp_sep), "List of Functional Suppression Profiles per Environment")
  
  # Test summary method for the list
  summ_sep <- summary(fsp_sep)
  expect_s3_class(summ_sep, "summary.functional_suppression_profiles_list")
  expect_length(summ_sep, 2)
  expect_output(print(summ_sep), "Environment: Env1")
  
  # Test plot method for the list
  plots_sep <- plot(fsp_sep, type = "profiles")
  expect_type(plots_sep, "list")
  expect_length(plots_sep, 2)
  expect_s3_class(plots_sep$Env1, "ggplot")
  # Test Joint Mode with response_scale = "percent" to ensure ... works and no warning is thrown
  expect_no_warning({
    functional_suppression_profiles(
      data = sim_dat_multi,
      reference = "Control",
      time = "time",
      response = "severity",
      treatment = "treatment",
      environment = "env",
      by_environment = FALSE,
      threshold = 5,
      k = 2,
      response_scale = "percent"
    )
  })
})
