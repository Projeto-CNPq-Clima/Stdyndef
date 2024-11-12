#' Title
#'
#' @param output a value
#' @param response A matrix of observed values for the response variable, with dimensions `Txn`, where `T` is the length of the time series and `n` is the number of locations.
#' @param FT A covariate matrix with dimensions `pxn`, where `p` is the number of coefficients and `n` is the number of locations for time `t = 1`.
#' @param MatFFT A covariate matrix with dimensions `Tx(pxn)`. Each row contains the vectorized elements of the `FT` matrix for each time `t`.
#' @param GT An evolution matrix with dimensions `pxp`.
#' @param sites A matrix `nx2` containing the geographic coordinates of the locations.
#' @param GAMA A matrix with dimensions `qxr` such that `q * r = T`. Each row contains the positions of the set `AT = {1, ..., T}`, defining the subset `GAMA_i = AT[GAMA[i,]]` for `i = 1, ..., q`. The `GAMA` matrix must guarantee a partition of `AT`, meaning that `GAMA_i ∩ GAMA_j` is empty for all `i ≠ j`, and the union of all `GAMA_i` equals `AT`.
#'
#'
#'#' @param GAMA A matrix that defines subsets of time indices to partition the time series into `q` disjoint groups.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{DIC}{The Deviance Information Criterion, a measure of model fit and complexity. Lower values indicate a better balance between fit and complexity.}
#'   \item{Pd}{The effective number of parameters, representing the complexity of the model.}
#'   \item{mean_deviance}{The posterior mean deviance, reflecting the fit of the model to the observed data.}
#' }
#'
#' @details
#' The Deviance Information Criterion (DIC) is calculated based on the posterior samples obtained from
#' the `SpatialDeformationMCMC` function. It balances model fit and complexity, providing a measure
#' that penalizes models with higher complexity (more parameters) while favoring good fit to the data.
#' The DIC is particularly useful in a Bayesian framework, allowing for model comparison and selection.
#'
#' @examples
#' # Load required data and libraries
#' library(Stdyndef)
#' data(temperature)
#' data(FT)
#' data(MatFFT)
#' data(GT)
#' data(sites)
#' data(GAMA)
#'
#' # Run MCMC to obtain posterior estimates
#' Mod <- SpatialDeformationMCMC(
#'   response = temperature,
#'   FT = FT,
#'   MatFFT = MatFFT,
#'   GT = GT,
#'   sites = sites,
#'   GAMA = GAMA,
#'   iteration = 100,
#'   burnin = 50,
#'   jump = 1
#' )
#'
#' # Calculate DIC for the model
#' dic_results <- DIC(
#'   output = Mod,
#'   response = temperature,
#'   FT = FT,
#'   MatFFT = MatFFT,
#'   GT = GT,
#'   sites = sites,
#'   GAMA = GAMA
#' )
#'
#' # View results
#' dic_results$DIC          # Deviance Information Criterion
#' dic_results$Pd           # Effective number of parameters
#' dic_results$mean_deviance # Posterior mean deviance
#'
#' @export
DIC <- function(output, response, FT, MatFFT, GT, sites, GAMA) {
  # Extract parameters from the MCMC output
  MTheta <- output$MTheta
  Mkappa <- output$Mkappa
  Msigmak <- output$Msigmak
  MPhi <- output$MPhi
  MDef <- output$MDef

  # Estimate parameters
  Thetaest <- matrix(apply(MTheta, 2, mean), nrow(response), ncol(MTheta) / nrow(response))
  qq <- length(GAMA)  # Number of regions
  ind1 <- seq(1, 2 * ncol(response) * qq, 2 * ncol(response))
  ind2 <- seq(2 * ncol(response), 2 * ncol(response) * qq, 2 * ncol(response))

  DEFest <- apply(MDef, 2, mean)
  MDD <- NULL
  for (l in 1:qq) {
    MDD <- rbind(MDD, t(as.matrix(as.vector(matrix(DEFest[ind1[l]:ind2[l]], ncol(response), 2)))))
  }

  kaest <- apply(Mkappa, 2, mean)
  vest <- apply(Msigmak, 2, mean)
  best <- apply(MPhi, 2, mean)

  # Calculate the likelihood for DIC
  DDpost <- -2 * logvero(MatFFT, response, Thetaest, kaest, vest, best, MDD, GAMA)

  # Calculate pD
  diffq <- NULL
  for (i in 1:nrow(Mkappa)) {
    kai <- Mkappa[i, ]
    vei <- Msigmak[i, ]
    bei <- MPhi[i, ]
    Thetaj <- matrix(MTheta[i, ], nrow(response), ncol(MTheta) / nrow(response))

    DEFj <- MDef[i, ]
    MDDj <- NULL
    for (l in 1:qq) {
      MDDj <- rbind(MDDj, t(as.matrix(as.vector(matrix(DEFj[ind1[l]:ind2[l]], ncol(response), 2)))))
    }

    temp1 <- -2 * logvero(MatFFT, response, Thetaj, kai, vei, bei, MDDj, GAMA)
    diffq <- c(diffq, temp1)
  }

  pD <- mean(diffq) - DDpost
  DIC_value <- mean(diffq) + pD

  return(list(DIC = DIC_value, pD = pD, DDpost = DDpost))
}
