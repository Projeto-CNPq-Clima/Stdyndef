#' Spatial Interpolation based on MCMC Spatiotemporal Model Output
#'
#' This function performs spatial interpolation for non-monitored locations based on the output of the `SpatialDeformationMCMC` model.
#'
#' @param MCMC_output A list containing the output of `SpatialDeformationMCMC`, including `MDef`, `MTheta`, `Msigmak`, `Mkappa`, and `MPhi`.
#' @param response A matrix of observed values for the response variable.
#' @param FT A matrix of covariates for the observed locations.
#' @param MatFFT A covariate matrix across time and space.
#' @param GT An evolution matrix.
#' @param sites A matrix containing the coordinates of observed locations.
#' @param sitesNO A matrix containing the coordinates of non-monitored locations to interpolate.
#' @param indi A vector of indices representing subsets of sites.
#' @param FTNO A matrix of covariates for the non-monitored locations.
#' @param MAP A matrix representing the map area.
#' @param bd A scalar parameter for the spatial correlation function.
#' @param dd Integer indicating the grid resolution (default: 11).
#' @return A matrix, `MatvecNO`, containing the interpolated values at the non-monitored locations for each time step.
#'
#'@export

spatial_interpolation <- function(MCMC_output, response, FT, MatFFT, GT, sites, sitesNO, indi, FTNO, MAP, bd, dd = 11) {

  # Extract elements from MCMC output
  MDef <- MCMC_output$MDef
  MTheta <- MCMC_output$MTheta
  Msigmak <- MCMC_output$Msigmak
  Mkappa <- MCMC_output$Mkappa
  MPhi <- MCMC_output$MPhi

  # Define grid for the area
  minx <- min(MAP[, 1])
  maxx <- max(MAP[, 1])
  miny <- min(MAP[, 2])
  maxy <- max(MAP[, 2])
  se1 <- seq((minx + (maxx - minx) / (2 * dd)), (maxx - (maxx - minx) / (2 * dd)), length.out = dd)
  se2 <- seq((miny + (maxy - miny) / (2 * dd)), (maxy - (maxy - miny) / (2 * dd)), length.out = dd)

  # Create grid matrix
  MTCELL <- NULL
  for (i in 1:length(se1)) {
    MCELL <- matrix(NA, nrow = length(se2), ncol = 2)
    for (j in 1:length(se2)) {
      MCELL[j, ] <- c(se1[i], se2[j])
    }
    MTCELL <- rbind(MTCELL, MCELL)
  }

  # Combine grid with observed sites
  DTotal <- rbind(MTCELL, sites)
  STotal <- complex(real = DTotal[, 1], imaginary = DTotal[, 2])
  stemp1 <- H(nrow(DTotal)) %*% STotal
  SzT <- t(H(nrow(DTotal))) %*% stemp1
  SzT <- cbind(Re(SzT), Im(SzT))

  Sz1 <- SzT[(nrow(MTCELL) + 1):nrow(SzT), ]
  Sz2 <- SzT[1:nrow(MTCELL), ][indi, ]
  D <- Sz1 - SzT[(nrow(MTCELL) + 1):nrow(SzT), ]
  Sz2c <- cbind(rep(D[1, 1], nrow(Sz2)), rep(D[1, 2], nrow(Sz2))) + Sz2

  VecSNO <- as.vector(Sz2c)
  Szz1 <- rbind(Sz2c, Sz1)
  distancias <- as.matrix(dist(Szz1))
  Rd1 <- exp(-bd * distancias)
  RA12 <- Rd1[1:nrow(Sz2c), (nrow(Sz2c) + 1):ncol(Rd1)]

  # Initialize result matrix
  MatvecNO <- NULL
  vecm <- as.vector(Sz1)
  vecmu <- as.vector(Sz2c)

  # Loop through MCMC samples
  for (l in 1:nrow(MTheta)) {
    DD <- NULL
    for (jj in 1:(ncol(MDef) / length(sites))) {
      ind1 <- (2 * ncol(response) * (jj - 1) + 1)
      ind2 <- (2 * ncol(response) * jj)
      DD <- cbind(DD, MDef[l, ind1:ind2])
    }

    VecDNO <- NULL
    mm <- rbind(Sz1, Sz2c)
    distancias <- as.matrix(dist(mm))
    RD <- exp(-bd * distancias)
    Rd <- RD[1:ncol(response), 1:ncol(response)]
    Rdug <- RD[(ncol(response) + 1):nrow(mm), 1:ncol(response)]
    I2 <- diag(2)

    for (jj in 1:(ncol(MDef) / length(sites))) {
      rtem <- vecmu + kronecker(I2, Rdug %*% solve(Rd)) %*% (DD[, jj] - vecm)
      VecDNO <- cbind(VecDNO, rtem)
    }

    # Perform interpolation for each time step
    Betaj <- matrix(MTheta[l, ], ncol = ncol(FT), byrow = TRUE)
    vecyNO <- NULL

    for (t in 1:nrow(response)) {
      Dzl <- rbind(matrix(VecDNO[, indi[t]], ncol = 2), matrix(DD[, indi[t]], ncol = 2))
      SIGMAl <- Msigmak[l, indi[t]] * ((1 - Mkappa[l, indi[t]]) * gCorr(MPhi[l, indi[t]], Dzl) + Mkappa[l, indi[t]] * diag(nrow(Dzl)))
      SIGM <- SIGMAl[(nrow(Sz2c) + 1):ncol(SIGMAl), (nrow(Sz2c) + 1):ncol(SIGMAl)]
      SIGMNO <- SIGMAl[1:nrow(Sz2c), 1:nrow(Sz2c)]
      SIGEst <- t(SIGMAl[(nrow(Sz2c) + 1):ncol(SIGMAl), 1:nrow(Sz2c)])

      mut <- t(FT) %*% Betaj[t, ]
      muNOt <- FTNO %*% Betaj[t, ]
      A <- muNOt + SIGEst %*% solve(SIGM) %*% (response[t, ] - mut)
      B <- SIGMNO - SIGEst %*% solve(SIGM) %*% t(SIGEst)

      yNO <- mvrnorm(1, A, B)
      vecyNO <- cbind(vecyNO, yNO)
    }

    MatvecNO <- rbind(MatvecNO, vecyNO)
  }

  return(MatvecNO)
}
