#' Border Map of Southern Brazil
#'
#' The `Map` dataset provides the geographic coordinates outlining the border of the Southern region of Brazil.
#' This dataset is useful for spatial visualizations and mapping purposes, especially in analyses focusing on
#' this specific geographic area.
#'
#' @format A data frame with `m` rows and 2 variables:
#' \describe{
#'   \item{longitude}{A numeric vector representing the longitude of each point on the border outline, in decimal degrees.}
#'   \item{latitude}{A numeric vector representing the latitude of each point on the border outline, in decimal degrees.}
#' }
#'
#' @details
#' The `Map` dataset captures the geographic boundary of the Southern region of Brazil, defined by a series
#' of latitude and longitude coordinates that trace the contour of the region. This data can be used as a
#' background layer for mapping applications or as a reference for spatial analyses within this region.
#'
#' @source \url{https://portal.inmet.gov.br/} (National Institute of Meteorology, Brazil)
#'
#' @examples
#' # Load the Map dataset
#' data(Map)
#'
#' # Plot the map outline, assuming Map is a matrix
#' plot(Map[,1], Map[,2], type = "l",
#'      main = "Border Map of Southern Brazil",
#'      xlab = "Longitude", ylab = "Latitude")
#'
"Map"
