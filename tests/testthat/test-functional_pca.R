test_that("functional_pca rejects invalid object classes", {
  invalid_obj <- list(curves = data.frame(time = 1:5, mu = 1:5))
  class(invalid_obj) <- "some_other_class"
  expect_error(functional_pca(invalid_obj), "must be of class 'functional_curves'")
})

test_that("functional_pca errors clearly when fitted curves are unavailable", {
  invalid_fc <- list(vars = list(time = "time", treatment = "trt"))
  class(invalid_fc) <- "functional_curves"
  expect_error(functional_pca(invalid_fc), "Fitted curves are not available")
})

test_that("functional_pca works properly on valid input", {
  # Mock a functional_curves object
  t_grid <- seq(1, 10, length.out = 10)
  trts <- c("A", "B", "C")
  
  curves_df <- expand.grid(time = t_grid, trt = trts)
  # Simple curves
  curves_df$mu <- ifelse(curves_df$trt == "A", curves_df$time * 0.1, 
                         ifelse(curves_df$trt == "B", curves_df$time * 0.2, curves_df$time * 0.3))
  
  mock_fc <- list(
    curves = curves_df,
    grid = t_grid,
    vars = list(time = "time", treatment = "trt")
  )
  class(mock_fc) <- "functional_curves"
  
  res <- functional_pca(mock_fc, center = TRUE, scale = FALSE)
  
  # 3. output inherits from "r4pde_functional_pca"
  expect_s3_class(res, "r4pde_functional_pca")
  
  # 4. scores have one row per curve
  expect_equal(nrow(res$scores), length(trts))
  expect_true("curve_id" %in% names(res$scores))
  
  # 5. eigenfunctions have values over the full time grid
  ef <- res$eigenfunctions
  expect_equal(length(unique(ef$time)), length(t_grid))
  expect_true(all(t_grid %in% ef$time))
  
  # 6. variance proportions sum to approximately 1 across all components
  # If we set n_components = NULL, var_explained = 1.0, it will return all.
  res_all <- functional_pca(mock_fc, center = TRUE, scale = FALSE, var_explained = 1.0)
  expect_equal(sum(res_all$variance$prop_var), 1, tolerance = 1e-5)
  
  # 7. retained components satisfy n_components or var_explained
  res_n1 <- functional_pca(mock_fc, n_components = 1)
  expect_equal(nrow(res_n1$variance), 1)
  expect_equal(res_n1$settings$n_retained, 1)
  
  # 8. reconstruction works and returns the expected columns
  expect_true(all(c("curve_id", "time", "fitted", "reconstructed") %in% names(res$reconstructed)))
  
  # 9. extractor functions return data frames
  scores_df <- get_fpca_scores(res)
  expect_s3_class(scores_df, "data.frame")
  
  ef_df <- get_fpca_eigenfunctions(res)
  expect_s3_class(ef_df, "data.frame")
  
  var_df <- get_fpca_variance(res)
  expect_s3_class(var_df, "data.frame")
  
  aug_df <- augment_functional_pca(res)
  expect_s3_class(aug_df, "data.frame")
  expect_true("residual" %in% names(aug_df))
  
  recon_df <- reconstruct_curves(res, components = 1)
  expect_s3_class(recon_df, "data.frame")
  
  # 10. plot methods return ggplot objects if ggplot2 is used
  p1 <- plot(res, type = "scree")
  expect_s3_class(p1, "ggplot")
  
  p2 <- plot(res, type = "components")
  expect_s3_class(p2, "ggplot")
  
  p3 <- plot(res, type = "reconstruction")
  expect_s3_class(p3, "ggplot")
})
