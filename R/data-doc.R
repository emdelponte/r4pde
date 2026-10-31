#' BlastWheat dataset
#'
#' Wheat blast dataset with severity and weather covariates.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{heading}{Date of heading}
#'   \item{inc_mean}{Mean incidence}
#'   \item{index_mean}{Mean wheat blast index}
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
#'   \item{i}{Number of diseased plants in the quadrat (used as \code{i} by \code{BPL()})}
#'   \item{n}{Number of plants assessed in the quadrat (used as \code{n} by \code{BPL()})}
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

#' FHBWheatBrazil dataset
#'
#' Fusarium head blight (FHB) index in wheat in Brazil.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{study}{Study ID or code}
#'   \item{year}{Crop year}
#'   \item{location}{Experimental site name}
#'   \item{state}{Brazilian state}
#'   \item{lat}{Latitude coordinate}
#'   \item{lon}{Longitude coordinate}
#'   \item{cultivar}{Wheat cultivar}
#'   \item{planting_date}{Planting date}
#'   \item{flowering}{Flowering date}
#'   \item{index}{FHB index}
#' }
#' @source Del Ponte Lab internal data
"FHBWheatBrazil"
#' SADRustSoybean dataset
#'
#' Severity of soybean rust estimated by 40 raters on 50 leaves, first without
#' and then with a standard area diagram (SAD). Raters 1 to 20 used the old SAD
#' and raters 21 to 40 the new SAD proposed by Franceschi et al. (2020). Use it
#' with [sad_data()] and the other `sad_*` functions.
#'
#' @format A data frame with 2000 rows and 6 columns:
#' \describe{
#'   \item{sad}{SAD used by the rater in the aided assessment: \code{"Old"} or \code{"New"}}
#'   \item{rater}{Rater identifier (1 to 40)}
#'   \item{leaf}{Leaf identifier (1 to 50)}
#'   \item{actual}{Actual severity (%) measured by image analysis}
#'   \item{unaided}{Severity (%) estimated without the SAD}
#'   \item{aided}{Severity (%) estimated with the SAD}
#' }
#' @source Franceschi, V. T., Alves, K. S., Mazaro, S. M., Godoy, C. V.,
#'   Duarte, H. S. S., and Del Ponte, E. M. (2020). A new standard area diagram
#'   set for assessment of severity of soybean rust improves accuracy of
#'   estimates and optimizes resource use. Plant Pathology 69:495-505.
"SADRustSoybean"
