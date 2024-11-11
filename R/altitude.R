#' Altitude of Monitoring Sites in Southern Brazil
#'
#' The `altitude` variable provides the altitude (in meters) for each of the 30 meteorological monitoring sites
#' located in the Southern region of Brazil. This data complements the geographic coordinates provided in the
#' `sites` dataset, offering additional spatial context for each location.
#'
#' @format A numeric vector with 30 elements:
#' \describe{
#'   \item{altitude}{Altitude of each monitoring site in meters above sea level.}
#' }s
#'
#' @details
#' The `altitude` variable indicates the elevation of each of the 30 monitoring stations used for collecting
#' meteorological data in Southern Brazil. This information can be used in conjunction with the temperature
#' data to assess the influence of altitude on climate variables across different locations.
#'
#' @source \url{https://portal.inmet.gov.br/}
#'
#' @examples
#' # Load the Altitude data
#' data(altitude)
#'
#' # Summary of altitude data
#' summary(altitude)
#'
#' # Plot altitude distribution
#' hist(altitude, main = "Distribution of Altitudes for Monitoring Sites",
#'      xlab = "Altitude (meters)", ylab = "Frequency")
#'
"altitude"
