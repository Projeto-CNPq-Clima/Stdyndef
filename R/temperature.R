#' Monthly Average Temperature Data for Southern Brazil (2007-2021)
#'
#' The `temperature` dataset provides the monthly average temperatures from January 2007 to December 2021
#' across 30 locations in the Southern region of Brazil. This dataset was compiled using publicly available
#' data from the National Institute of Meteorology (INMET) in Brazil.
#'
#' @format A data frame with 180 rows and 30 variables:
#' \describe{
#'   \item{rows}{Each row represents a specific month, from January 2007 to December 2021, covering a total of 180 months.}
#'   \item{columns}{Each column represents the monthly average temperature (in degrees Celsius) for one of the 30 locations in Southern Brazil.}
#' }
#'
#' @details
#' The dataset was constructed using historical data from INMET (the National Institute of Meteorology, Brazil).
#' It includes monthly average temperature data for 30 monitoring stations located across the Southern region
#' of Brazil. The data spans from January 2007 to December 2021, with each column containing the time series
#' for a specific location.
#'
#' @source \url{https://portal.inmet.gov.br/}
#'
#' @examples
#' # Load the dataset
#' data(temperature)
#'
#' # View the first few rows of the dataset
#' head(temperature)
#'
#' # Plot the average temperature for the first location
#' plot(temperature[, 1], type = "l", xlab = "Month", ylab = "Temperature (°C)",
#'      main = "Monthly Average Temperature for Location 1")
#'
"temperature"
