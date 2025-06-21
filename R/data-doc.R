#' BlastWheat dataset
#'
#' Wheat blast dataset with severity and weather covariates.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{heading}{Date of heading}
#'   \item{inc_mean}{Mean incidence}
#'   \item{index_mean}{FHB index mean}
#'   \item{latitude}{Latitude coordinate}
#'   \item{location}{Experimental site name}
#'   \item{longitude}{Longitude coordinate}
#'   \item{state}{Brazilian state}
#'   \item{study}{Study ID or code}
#'   \item{year}{Crop year}
#'   \item{yld_mean}{Mean yield}
#' }
#' @source Del Ponte Lab internal data
"BlastWheat"
#' BudBlightSoybean dataset
#'
#' Soybean bud blight incidence in experimental blocks.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{block}{Block number}
#'   \item{time}{Time point of assessment}
#'   \item{treat}{Treatment name}
#'   \item{y}{Incidence or severity value}
#' }
#' @source Del Ponte Lab internal data
"BudBlightSoybean"
#' DidymellaWatermelon dataset
#'
#' Assessment of Didymella symptoms in watermelon plots.
#'
#' @format A data frame with:
#' \describe{
#'   \item{EW_row}{Row position (east–west)}
#'   \item{NS_col}{Column position (north–south)}
#'   \item{dap}{Days after planting}
#'   \item{severity}{Disease severity}
#' }
#' @source Del Ponte Lab internal data
"DidymellaWatermelon"
#' FHBWheat dataset
#'
#' Fusarium head blight quadrat assessments in wheat.
#'
#' @format A data frame with:
#' \describe{
#'   \item{field}{Field identifier}
#'   \item{i}{Row position}
#'   \item{n}{Column position}
#'   \item{quadrat}{Quadrat ID}
#'   \item{season}{Crop season}
#' }
#' @source Del Ponte Lab internal data
"FHBWheat"
#' FusariumBanana dataset
#'
#' Observations of Fusarium symptoms in banana fields.
#'
#' @format A data frame with:
#' \describe{
#'   \item{field}{Field ID}
#'   \item{lat}{Latitude}
#'   \item{lon}{Longitude}
#'   \item{marker}{Infection marker presence}
#' }
#' @source Del Ponte Lab internal data
"FusariumBanana"
#' RustSoybean dataset
#'
#' Soybean rust severity and field metadata.
#'
#' @format A data frame with:
#' \describe{
#'   \item{detection}{Detection score or date}
#'   \item{epidemia}{Epidemic phase}
#'   \item{latitude}{Latitude}
#'   \item{local}{Location name}
#'   \item{longitude}{Longitude}
#'   \item{planting}{Planting date or stage}
#'   \item{severity}{Disease severity}
#' }
#' @source Del Ponte Lab internal data
"RustSoybean"
#' WhiteMoldSoybean dataset
#'
#' National dataset of white mold severity and yield.
#'
#' @format A data frame with:
#' \describe{
#'   \item{country}{Country name}
#'   \item{elevation}{Field elevation}
#'   \item{elevation_class}{Elevation class}
#'   \item{harvest_year}{Year of harvest}
#'   \item{inc}{Incidence}
#'   \item{inc_check}{Check plot incidence}
#'   \item{inc_class}{Incidence class}
#'   \item{location}{Location name}
#'   \item{region}{Geographical region}
#'   \item{scl}{Soybean canopy layer}
#'   \item{season}{Crop season}
#'   \item{state}{State name}
#'   \item{study}{Study identifier}
#'   \item{treat}{Treatment applied}
#'   \item{yld}{Yield}
#'   \item{yld_check}{Yield of untreated check}
#'   \item{yld_class}{Yield class}
#' }
#' @source Del Ponte Lab internal data
"WhiteMoldSoybean"
#' SpatialAggregated dataset
#'
#' Simulated aggregated spatial binary disease pattern.
#'
#' @format A data frame with:
#' \describe{
#'   \item{x}{x-coordinate}
#'   \item{y}{y-coordinate}
#' }
#' @source Simulated example
"SpatialAggregated"

#' SpatialRandom dataset
#'
#' Simulated random spatial binary disease pattern.
#'
#' @format A data frame with:
#' \describe{
#'   \item{x}{x-coordinate}
#'   \item{y}{y-coordinate}
#' }
#' @source Simulated example
"SpatialRandom"


