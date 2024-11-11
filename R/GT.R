#' Transition Matrix for State-Space Model at Time t
#'
#' The `GT` matrix is the transition matrix used in the state-space model at a given time `t`. This matrix
#' describes the evolution of the state variables over time, allowing the model to capture the temporal
#' dependencies and transitions between states in the spatiotemporal process.
#'
#' @format A matrix with `p` rows and `p` columns:
#' \describe{
#'   \item{p}{The number of state variables in the model.}
#' }
#'
#' @details
#' The `GT` matrix serves as the transition matrix in the state-space model at time `t`. It defines how the
#' state variables at time `t - 1` influence the state variables at time `t`, playing a crucial role in
#' modeling the temporal dynamics of the system. This matrix is essential for capturing the evolution and
#' propagation of states across time steps in the model.
#'
#' @examples
#' # Load the GT matrix
#' data(GT)
#'
#' # Display the structure of GT
#' str(GT)
#'
"GT"
