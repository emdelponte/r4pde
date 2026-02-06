#' Fetch BR-DWGD Data for Multiple Locations
#'
#' This function extracts daily weather data from the Brazilian Daily Weather Gridded Data (BR-DWGD)
#' NetCDF files for specified locations and dates. It mimics the functionality of `get_nasapower` 
#' and `get_era5` by returning a daily-aggregated data frame with columns named according 
#' to `nasapower` conventions.
#'
#' @param data A data frame containing the input data, including columns for `latitude`, `longitude`, 
#' and the date column specified by `date_col`.
#' @param days_around An integer specifying the number of days before and after the date in the date column
#' to download data.
#' @param date_col A character string specifying the name of the date column in the data frame.
#' @param study_col A character string specifying the name of the column containing the study identifier 
#' in the input data frame (default: "study").
#' @param path A character string specifying the directory where the BR-DWGD NetCDF files are stored.
#' @param vars A character vector specifying the weather variables to fetch. 
#' These are expected to be the prefix of the NetCDF files (e.g., "pr", "Tmax").
#' Default: `c("pr", "Tmax", "Tmin", "Rs", "RH")`.
#'
#' @return A data frame (tibble) with daily weather data. Columns are named to match 
#' `nasapower` conventions where possible: `date`, `PRECTOTCORR` (from `pr`), 
#' `T2M_MAX` (from `Tmax`), `T2M_MIN` (from `Tmin`), `ALLSKY_SFC_SW_DWN` (from `Rs`), 
#' `RH2M` (from `RH`), `longitude`, `latitude`, and `study`.
#'
#' @details
#' The function requires the `terra` package to read NetCDF files and the `tidyr` package 
#' for data reshaping. It expects NetCDF files to be named starting with the variable name 
#' (e.g., `pr_*.nc`) and to contain the variable with the same name.
#' 
#' @importFrom dplyr mutate bind_rows as_tibble select everything inner_join
#' @importFrom lubridate days as_date
#' @importFrom purrr map_dfr
#' @importFrom progress progress_bar
#' @importFrom tidyr pivot_wider
#' @importFrom terra rast time extract cellFromXY nlyr
#' @importFrom rlang .data
#' @export
#' @family Disease modeling
get_brdwgd <- function(data, days_around, date_col, study_col = "study",
                      path = "netcdf_files", 
                      vars = c("pr", "Tmax", "Tmin", "Rs", "RH")) {
  
  # 1. Flexible coordinate and study column detection
  lon_col <- if ("longitude" %in% names(data)) "longitude" else if ("lon" %in% names(data)) "lon" else stop("Longitude column not found.")
  lat_col <- if ("latitude" %in% names(data)) "latitude" else if ("lat" %in% names(data)) "lat" else stop("Latitude column not found.")
  
  # 2. Pre-calculate all target dates for all locations to find the global range
  data$..row_id <- 1:nrow(data)
  date_list <- lapply(1:nrow(data), function(i) {
    base_date <- as.Date(data[[date_col]][i])
    seq(base_date - days_around, base_date + days_around, by = "day")
  })
  all_needed_dates <- unique(do.call(c, date_list))
  
  # 3. Load and cache search for variables
  res_all_vars <- list()
  
  # Mapping of BR-DWGD variable names to NASA POWER/ERA5 conventions
  var_name_map <- c(
    pr = "PRECTOTCORR",
    Tmax = "T2M_MAX",
    Tmin = "T2M_MIN",
    Rs = "ALLSKY_SFC_SW_DWN",
    RH = "RH2M"
  )

  for (v in vars) {
    pattern <- paste0("^", v, ".*\\.nc$")
    files <- list.files(path, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) {
      warning(paste("No files for", v))
      next
    }
    
    # Message for transparency
    message(paste("Processing", v, "..."))
    
    # Load once
    r <- terra::rast(files)
    t_vals <- terra::time(r)
    if (is.numeric(t_vals)) t_vals <- as.Date(t_vals, origin = "1961-01-01") else t_vals <- as.Date(t_vals)
    
    # Identify unique layers needed for this variable across ALL locations
    global_idx <- which(t_vals %in% all_needed_dates)
    if (length(global_idx) == 0) next
    
    # Subset once
    r_sub <- r[[global_idx]]
    t_sub <- t_vals[global_idx]
    
    # Extract ALL points at once (Matrix extraction is much faster than point by point)
    pts <- matrix(c(data[[lon_col]], data[[lat_col]]), ncol = 2)
    # This returns a matrix: rows = locations, cols = layers
    extracted_matrix <- terra::extract(r_sub, pts)
    
    # terra::extract might return a data frame with an ID column if pts is a data frame or SpatVector,
    # but for a matrix it usually returns a matrix. Let's ensure we have just the values.
    if (inherits(extracted_matrix, "data.frame") && "ID" %in% names(extracted_matrix)) {
      extracted_matrix <- as.matrix(extracted_matrix[, -1, drop = FALSE])
    } else {
      extracted_matrix <- as.matrix(extracted_matrix)
    }
    
    # Map variable name
    mapped_name <- if (v %in% names(var_name_map)) var_name_map[[v]] else v
    
    # Organize back to a long data frame
    var_results <- list()
    for (i in 1:nrow(data)) {
      # Target dates for this specific row
      row_dates <- date_list[[i]]
      # Column indices in the extracted_matrix that correspond to these dates
      col_idx <- which(t_sub %in% row_dates)
      
      if (length(col_idx) > 0) {
        var_results[[i]] <- data.frame(
          ..row_id = i,
          date = t_sub[col_idx],
          value = as.numeric(extracted_matrix[i, col_idx]),
          var = mapped_name,
          stringsAsFactors = FALSE
        )
      }
    }
    res_all_vars[[v]] <- do.call(rbind, var_results)
  }

  if (length(res_all_vars) == 0) return(dplyr::as_tibble(data.frame()))

  # Combine all variables
  combined_long <- do.call(rbind, res_all_vars)
  
  # Final Join with original data to get study and coords
  final_df <- combined_long %>%
    dplyr::inner_join(data, by = "..row_id") %>%
    tidyr::pivot_wider(names_from = "var", values_from = "value") %>%
    dplyr::select(-"..row_id") %>%
    dplyr::mutate(study = .data[[study_col]])

  # Cleanup selection
  final_res <- final_df %>%
    dplyr::select(dplyr::any_of(c("date", "study", "latitude", "longitude", lon_col, lat_col)), 
                  dplyr::everything())

  return(dplyr::as_tibble(final_res))
}

