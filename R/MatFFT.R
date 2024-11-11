#' Covariate Matrix for State-Space Model at Each Time Step (Vectorized)
#'
#' The `MatFFT` matrix provides the covariates for the state-space model across multiple time steps,
#' where each row represents the `FT` matrix (covariate matrix at a specific time `t`) in a vectorized form.
#' This structure enables efficient handling of time-dependent covariates across all observed time steps.
#'
#' @format A matrix with `T` rows and `p * n` columns:
#' \describe{
#'   \item{T}{The number of time steps in the model.}
#'   \item{p * n}{The number of covariate elements for each time step, where `p` is the number of covariates
#'   and `n` is the number of spatial locations. Each row is a vectorized representation of the `FT` matrix
#'   at a given time step.}
#' }
#'
#' @details
#' The `MatFFT` matrix is constructed by vectorizing the `FT` matrix for each time step and stacking these
#' vectorized matrices row-wise. Thus, the first row contains the vectorized `FT` matrix for time `t = 1`,
#' the second row contains the vectorized `FT` matrix for time `t = 2`, and so on. This structure facilitates
#' the incorporation of time-varying covariates in the state-space model.
#'
#' @examples
#' # Load the MatFFT matrix
#' data(MatFFT)
#'
#' # Display the structure of MatFFT
#' str(MatFFT)
#'
"MatFFT"
