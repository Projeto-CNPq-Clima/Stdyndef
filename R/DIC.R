#' Deviance Information Criterion (DIC) for Spatiotemporal Model with Dynamic Deformation
#'
#' Calculates the Deviance Information Criterion (DIC) for a spatiotemporal model with dynamic deformation
#' and nonstationary covariance structures. The DIC is used to assess model fit while penalizing model
#' complexity, facilitating model comparison and selection.
#'
#' @param MCMC_output A list containing the output from the `SpatialDeformationMCMC` function, which includes
#' the posterior samples needed to compute the deviance.
#' @param response A matrix of observed values for the response variable across multiple time points and locations.
#' @param FT A covariate matrix for the model at time `t = 1`.
#' @param MatFFT A matrix containing vectorized covariate matrices for each time step, where each row represents
#' the covariates at a given time in a vectorized form.
#' @param GT An evolution matrix used in the state-space model.
#' @param sites A matrix containing the geographic coordinates of the observation locations.
#' @param GAMA A matrix that defines subsets of time indices to partition the time series into `q` disjoint groups.
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
#'   MCMC_output = Mod,
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
DIC <- function(MCMC_output, response, FT, MatFFT, GT, sites, GAMA) {
  # Extract parameters from the MCMC output
  MTheta <- MCMC_output$MTheta
  Mkappa <- MCMC_output$Mkappa
  Msigmak <- MCMC_output$Msigmak
  MPhi <- MCMC_output$MPhi
  MDef <- MCMC_output$MDef

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
