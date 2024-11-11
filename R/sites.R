#' Geographic Coordinates of Monitoring Sites in Southern Brazil
#'
#' The `sites` dataset provides the geographic coordinates of 30 locations in the Southern region of Brazil.
#' These locations correspond to the meteorological monitoring stations used to collect monthly average temperature data.
#' The dataset was compiled using data from the National Institute of Meteorology (INMET) in Brazil.
#'
#' @format A data frame with 30 rows and 2 variables:
#' \describe{
#'   \item{latitude}{A numeric vector representing the latitude of each location in decimal degrees.}
#'   \item{longitude}{A numeric vector representing the longitude of each location in decimal degrees.}
#' }
#'
#' @details
#' The `sites` dataset contains the geographic coordinates (latitude and longitude) of 30 monitoring stations
#' located in the states of Paraná, Santa Catarina, and Rio Grande do Sul in Southern Brazil. These coordinates
#' can be used for spatial analysis or to visualize the spatial distribution of the meteorological stations.
#'
#' @source \url{https://portal.inmet.gov.br/}
#'
#' @examples
#' # Load the dataset
#' data(sites)
#'
#' # View the first few rows of the dataset
#' head(sites)
#'
"sites"
