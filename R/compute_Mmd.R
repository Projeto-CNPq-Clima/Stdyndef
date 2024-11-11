#' Computes the Mmd matrix using Spatial Deformation MCMC output
#'
#' @param MCMC_output A list containing the output of `SpatialDeformationMCMC`, which includes the matrix `MDef`.
#' @param response A matrix of observed values for the response variable.
#' @param sites A matrix containing the geographic coordinates of the locations.
#' @param Map A matrix representing the geographic coordinates of the map boundaries.
#' @param dd An integer specifying the resolution of the grid, where the estimated deformations will be
#' computed for a grid of dimension `(dd+1) x (dd+1)`.
#'
#' @return A matrix where each column represents the posterior estimate of the q-th dynamic deformation
#' on a grid of dimension `(dd+1) x (dd+1)`. Additionally, it includes the `Map` outline as a reference
#' for spatial context.
#'
#' @details
#' This function calculates the dynamic deformations by using posterior samples obtained from MCMC estimation.
#' For each posterior sample, the deformations are interpolated over a specified grid resolution `(dd+1) x (dd+1)`,
#' with each column in the resulting matrix corresponding to the posterior deformation estimates at that grid
#' resolution. The `Map` outline is returned alongside these estimates for geographic reference.
#'
#' @examples
#' # Example usage of compute_Mmd
#' MCMC_output <- SpatialDeformationMCMC(response = temperature, FT = FT, MatFFT = MatFFT, GT = GT,
#'                                        sites = sites, GAMA = GAMA, iteration = 100, burnin = 50, jump = 1)
#' response <- temperature
#' sites <- sites
#' Map <- Map
#' dd <- 10
#'
#' Mmd_result <- compute_Mmd(MCMC_output, response, sites, Map, dd)
#'
#' # Plot the posterior deformations on the grid with Map outline
#' image(matrix(Mmd_result[,1], nrow = dd + 1), main = "Posterior Deformation (q = 1)")
#' lines(Map[,1], Map[,2], col = "red")
#'
#' @export
compute_Mmd <- function(MCMC_output, response, sites, Map, dd) {

  # Extract from MCMC output
  MDef <- MCMC_output$MDef
  bd<-MCMC_output$bd

  # Initialize parameters
  Mmd <- NULL
  qq <- ncol(MDef) / length(sites)
  n <- ncol(response)
  DEFest <- apply(MDef, 2, mean)
  ind1 <- seq(1, 2 * n * qq, 2 * n)
  ind2 <- seq(2 * n, 2 * n * qq, 2 * n)

  # Loop through each subset of deformations
  for (j in 1:qq) {
    D2est <- matrix(DEFest[ind1[j]:ind2[j]], n, 2)
    Da <- D2est

    # Combine map and site coordinates
    DTotal <- rbind(Map, sites)
    STotal <- complex(real = DTotal[, 1], imaginary = DTotal[, 2])
    stemp1 <- H(nrow(DTotal)) %*% STotal  # Remove translation
    SzT <- t(H(nrow(DTotal))) %*% stemp1  # Translated configuration
    SzT <- cbind(Re(SzT), Im(SzT))

    # Define grid
    minx <- min(SzT[, 1])
    maxx <- max(SzT[, 1])
    miny <- min(SzT[, 2])
    maxy <- max(SzT[, 2])
    se1 <- seq(minx, maxx, length.out = dd + 1)
    se2 <- seq(miny, maxy, length.out = dd + 1)
    MTCELL <- NULL

    # Create grid matrix
    for (i in 1:length(se1)) {
      MCELL <- matrix(NA, nrow = length(se2), ncol = 2)
      for (k in 1:length(se2)) {
        MCELL[k, ] <- c(se1[i], se2[k])
      }
      MTCELL <- rbind(MTCELL, MCELL)
    }

    # Create combined matrix for distance calculations
    matrizz <- rbind(MTCELL, SzT[1:nrow(Map), ])
    m <- SzT[(nrow(Map) + 1):nrow(SzT), ]
    mu <- matrizz
    n1 <- nrow(mu)

    vecm <- rbind(matrix(m[, 1], ncol = 1), matrix(m[, 2], ncol = 1))
    vecD <- rbind(matrix(Da[, 1], ncol = 1), matrix(Da[, 2], ncol = 1))
    vecmu <- rbind(matrix(mu[, 1], ncol = 1), matrix(mu[, 2], ncol = 1))

    mm <- rbind(m, mu)
    distancias <- as.matrix(dist(mm))
    RD <- exp(-bd * distancias)
    Rd <- RD[1:n, 1:n]
    Rdug <- RD[(n + 1):nrow(mm), 1:n]
    I2 <- diag(1, 2)

    # Calculate vecmd and update Mmd
    vecmd <- vecmu + (I2 %x% (Rdug %*% solve(Rd))) %*% (vecD - vecm)
    Mmd <- cbind(Mmd, as.matrix(vecmd))
  }

  return(Mmd)
}
