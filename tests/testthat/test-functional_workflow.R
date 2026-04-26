test_that("functional_curves returns expected class and components", {
  # Generate some dummy data
  set.seed(123)
  df <- data.frame(
    time = rep(1:5, 4),
    treatment = rep(c("A", "B"), each = 10),
    block = rep(rep(1:2, each = 5), 2),
    y = pmax(0, pmin(1, rep(c(0.1, 0.3, 0.5, 0.7, 0.9), 4) + rnorm(20, 0, 0.05)))
  )
  
  fc <- functional_curves(
    data = df, time = "time", response = "y", treatment = "treatment",
    block = "block", min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  expect_s3_class(fc, "functional_curves")
  expect_true(all(c("gam", "curves", "grid", "observed_data", "vars", "settings", "plot_mean", "scores") %in% names(fc)))
  expect_s3_class(fc$plot_mean, "ggplot")
})

test_that("functional_distances returns symmetric distance matrix with zero diagonal", {
  set.seed(123)
  df <- data.frame(
    time = rep(1:5, 4),
    treatment = rep(c("A", "B"), each = 10),
    block = rep(rep(1:2, each = 5), 2),
    y = pmax(0, pmin(1, rep(c(0.1, 0.3, 0.5, 0.7, 0.9), 4) + rnorm(20, 0, 0.05)))
  )
  fc <- functional_curves(
    data = df, time = "time", response = "y", treatment = "treatment",
    block = "block", min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  fd <- functional_distances(fc, cluster_k = 2, show_progress = FALSE)
  
  expect_s3_class(fd, "functional_distances")
  
  # Symmetric
  expect_equal(fd$distance_matrix, t(fd$distance_matrix))
  # Zero diagonal
  expect_true(all(diag(fd$distance_matrix) == 0))
})

test_that("functional_resistance rankings and instability penalties", {
  set.seed(123)
  # A is low disease, B is medium, C is high (Reference)
  df <- data.frame(
    time = rep(1:5, 6),
    treatment = rep(c("A", "B", "C"), each = 10),
    block = rep(rep(1:2, each = 5), 3),
    y = pmax(0, pmin(1, c(
      rep(c(0.01, 0.05, 0.1, 0.15, 0.2), 2), # A
      rep(c(0.1, 0.2, 0.3, 0.4, 0.5), 2),   # B
      rep(c(0.2, 0.4, 0.6, 0.8, 0.9), 2)    # C
    ) + rnorm(30, 0, 0.01)))
  )
  
  fc <- functional_curves(
    data = df, time = "time", response = "y", treatment = "treatment",
    block = "block", min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  fr_no_inst <- functional_resistance(fc, reference = "C", group_method = "quantile", n_groups = 3)
  expect_s3_class(fr_no_inst, "functional_resistance")
  
  # A should have lower disease than C, so higher FRI than B
  safri_a <- fr_no_inst$table$SAFRI[fr_no_inst$table$treatment == "A"]
  safri_b <- fr_no_inst$table$SAFRI[fr_no_inst$table$treatment == "B"]
  expect_true(safri_a > safri_b)
  
  # Instability penalty
  inst_mock <- data.frame(treatment = c("A", "B", "C"), nFI = c(1, 0, 0.5))
  fr_inst <- functional_resistance(fc, instability = inst_mock, reference = "C", lambda = 1)
  
  safri_a_inst <- fr_inst$table$SAFRI[fr_inst$table$treatment == "A"]
  # SAFRI should decrease when nFI increases (lambda=1)
  expect_true(safri_a_inst < safri_a)
})

test_that("compare_curves works as wrapper", {
  set.seed(123)
  df <- data.frame(
    time = rep(1:5, 4),
    treatment = rep(c("A", "B"), each = 10),
    block = rep(rep(1:2, each = 5), 2),
    y = pmax(0, pmin(1, rep(c(0.1, 0.3, 0.5, 0.7, 0.9), 4) + rnorm(20, 0, 0.05)))
  )
  suppressWarnings({
    cc <- compare_curves(
      data = df, time = "time", response = "y", treatment = "treatment",
      block = "block", min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
    )
  })
  expect_s3_class(cc, "r4pde_compare_curves")
  expect_true(!is.null(cc$gam))
  expect_true(!is.null(cc$distance))
})

test_that("functional_resistance bootstrap output", {
  set.seed(123)
  df <- data.frame(
    time = rep(1:5, 6),
    treatment = rep(c("A", "B", "C"), each = 10),
    block = rep(rep(1:2, each = 5), 3),
    y = pmax(0, pmin(1, c(
      rep(c(0.01, 0.05, 0.1, 0.15, 0.2), 2), # A
      rep(c(0.1, 0.2, 0.3, 0.4, 0.5), 2),   # B
      rep(c(0.2, 0.4, 0.6, 0.8, 0.9), 2)    # C
    ) + rnorm(30, 0, 0.01)))
  )
  fc <- functional_curves(
    data = df, time = "time", response = "y", treatment = "treatment",
    block = "block", min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  fr <- functional_resistance(fc, reference = "C", group_method = "bootstrap", n_boot = 50, n_groups = 3)
  
  # output returns one distribution per genotype
  expect_true(!is.null(fr$bootstrap$boot_distributions))
  dist_counts <- table(fr$bootstrap$boot_distributions$treatment)
  expect_equal(as.numeric(dist_counts["A"]), 50)
  
  # pairwise bootstrap differences separate clearly different genotypes
  # A and B are very different
  grp_a <- fr$table$bootstrap_group[fr$table$treatment == "A"]
  grp_b <- fr$table$bootstrap_group[fr$table$treatment == "B"]
  expect_true(grp_a != grp_b)
  
  # grouping labels are ordered from most resistant to least
  # A > B > C in resistance
  expect_true(as.numeric(factor(grp_a, levels = LETTERS)) < as.numeric(factor(grp_b, levels = LETTERS)))
})

test_that("functional_curves processes covariates correctly", {
  set.seed(123)
  df <- data.frame(
    time = rep(1:5, 4),
    treatment = rep(c("A", "B"), each = 10),
    block = rep(rep(1:2, each = 5), 2),
    covar = rep(c("early", "late"), each = 10),
    y = pmax(0, pmin(1, rep(c(0.1, 0.3, 0.5, 0.7, 0.9), 4) + rnorm(20, 0, 0.05)))
  )
  
  fc <- functional_curves(
    data = df, time = "time", response = "y", treatment = "treatment",
    block = "block", covariates = "covar", include_covariates = TRUE,
    covariate_smooths = TRUE, min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  expect_true(!is.null(fc$genotype_info))
  expect_true("covar" %in% names(fc$genotype_info))
  expect_equal(as.character(fc$genotype_info$covar), c("early", "late"))
  
  # Error if covariate varies within genotype
  df_bad <- df
  df_bad$covar[1] <- "late"
  expect_error(functional_curves(
    data = df_bad, time = "time", response = "y", treatment = "treatment",
    block = "block", covariates = "covar", include_covariates = TRUE,
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  ), "multiple values")
})

test_that("functional_resistance works with adjust_by and reference_within_group", {
  set.seed(123)
  df <- data.frame(
    time = rep(1:5, 8),
    treatment = rep(c("A_E", "B_E", "A_L", "B_L"), each = 10),
    block = rep(rep(1:2, each = 5), 4),
    covar = rep(c("E", "E", "L", "L"), each = 10),
    y = pmax(0, pmin(1, c(
      rep(c(0.1, 0.2, 0.3, 0.4, 0.5), 2), # A_E
      rep(c(0.2, 0.4, 0.6, 0.8, 0.9), 2), # B_E (Susc Early)
      rep(c(0.1, 0.2, 0.3, 0.4, 0.5), 2), # A_L
      rep(c(0.2, 0.4, 0.6, 0.8, 0.9), 2)  # B_L (Susc Late)
    ) + rnorm(40, 0, 0.01)))
  )
  fc <- functional_curves(
    data = df, time = "time", response = "y", treatment = "treatment",
    block = "block", covariates = "covar", include_covariates = FALSE,
    min_points = 3, show_progress = FALSE, family_try = "quasibinomial"
  )
  
  # adjust_by = covar, reference_within_group = TRUE
  refs <- c(E = "B_E", L = "B_L")
  fr <- functional_resistance(fc, reference = refs, adjust_by = "covar", reference_within_group = TRUE, group_method = "quantile", n_groups = 2, group_within_adjust_by = TRUE)
  
  expect_true(!is.null(fr$table$covar))
  expect_true("rank_within_group" %in% names(fr$table))
  
  # reference missing
  refs_bad <- c(E = "B_E", L = "Missing")
  expect_error(functional_resistance(fc, reference = refs_bad, adjust_by = "covar", reference_within_group = TRUE), "not found in group")
})
