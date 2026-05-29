#' Vietnam provinces (sf)
#'
#' Spatial polygons for Vietnamese provinces with administrative and climate
#' region codes.
#'
#' @format An `sf` data frame with 63 rows and 4 columns:
#' \describe{
#'   \item{region}{Administrative region code.}
#'   \item{province}{Province name.}
#'   \item{geometry}{Polygon geometry (sfc).}
#'   \item{climate_region}{Climate region code (may be `NA`).}
#' }
"vietnam_provinces"

#' Vietnam administrative regions
#'
#' Lookup table mapping administrative region codes to names.
#'
#' @format A data frame with 6 rows and 2 columns:
#' \describe{
#'   \item{code}{Region code.}
#'   \item{name}{Region name.}
#' }
"vietnam_regions"

#' Vietnam daily temperature records
#'
#' Annual records of daily minimum and maximum temperatures by province.
#'
#' @format A tibble with 1890 rows and 5 columns:
#' \describe{
#'   \item{year}{Year of observation.}
#'   \item{region}{Administrative region code.}
#'   \item{province}{Province name.}
#'   \item{t_max}{List-column of daily maximum temperatures.}
#'   \item{t_min}{List-column of daily minimum temperatures.}
#' }
"vietnam_temperature"

#' Vietnam temperatures as distributional data
#'
#' Same data as [vietnam_temperature] with `t_max` and `t_min` stored as
#' distributional (`ddl`) list-columns.
#'
#' @format A rowwise tibble with 1890 rows and 5 columns:
#' \describe{
#'   \item{year}{Year of observation.}
#'   \item{region}{Administrative region code.}
#'   \item{province}{Province name.}
#'   \item{t_max}{`ddl` list-column of daily maximum temperatures.}
#'   \item{t_min}{`ddl` list-column of daily minimum temperatures.}
#' }
"vietnam_temperature_dd"
