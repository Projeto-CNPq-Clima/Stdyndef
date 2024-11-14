#' Covariate Matrix for Unobserved Locations in the Spatiotemporal Model
#'
#' The `FTNO` matrix provides the covariates for locations where the spatiotemporal process was not observed.
#' This matrix is used during interpolation to estimate the process at these unobserved locations.
#'
#' @format A matrix with `p` rows and `n_unobs` columns:
#' \describe{
#'   \item{p}{The number of covariates or predictors in the model.}
#'   \item{n_unobs}{The number of spatial locations where the process was not observed.}
#' }
#'
#' @details
#' The `FTNO` matrix contains the covariate values corresponding to the locations where the spatiotemporal
#' process was not observed. These covariates are used in conjunction with the model to interpolate the
#' process at unobserved locations. Each row represents a covariate, and each column corresponds to a
#' specific unobserved location.
#'
#' This matrix is crucial for extending the model's predictions to areas without direct observations,
#' enabling spatial interpolation and prediction tasks.
#'
#' @examples
#' # Load the FTNO matrix
#' data(FTNO)
#'
#' # Display the structure of FTNO
#' str(FTNO)
#'
#' # Access the covariates for the first unobserved location
#' FTNO[, 1]
#'
#' # Plot the first covariate across all unobserved locations
#' plot(FTNO[1, ], type = "l", main = "First Covariate Across Unobserved Locations",
#'      xlab = "Unobserved Location Index", ylab = "Covariate Value")
#'
"FTNO"
