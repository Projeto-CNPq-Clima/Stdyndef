#' Spatio-Temporal Interpolation for Response Variable
#'
#' This function performs spatio-temporal interpolation for a response variable based on the output
#' of the SpatialDeformationMCMC model.
#'
#' @param MCMC_output A list containing the output matrices Msigmak, MPhi, Mkappa, MTheta, MDef, and the parameter bd from SpatialDeformationMCMC.
#' @param response Matrix of observed values for the response variable.
#' @param FT Covariate matrix for observed locations.
#' @param MatFFT Covariate matrix for all time-space points.
#' @param GT Evolution matrix.
#' @param sites Matrix of observed locations coordinates.
#' @param INDI Indices of the grid points for interpolation.
#' @param GAMA Matrix defining subsets of locations for interpolation.
#' @param FTNO Covariate matrix for non-observed locations.
#' @param RS Matrix representing the map area for spatial interpolation.
#' @return A matrix `MatvecNO` containing interpolated values for non-observed locations for each time step.
#' @export
spatio_temporal_interpolation <- function(MCMC_output, response, FT, MatFFT, GT, sites, INDI, GAMA, FTNO, RS) {

  # Extract matrices from MCMC output
  Msigmak <- MCMC_output$Msigmak
  MPhi <- MCMC_output$MPhi
  Mkappa <- MCMC_output$Mkappa
  MTheta <- MCMC_output$MTheta
  MDef <- MCMC_output$MDef
  bd <- MCMC_output$bd

  # Parameters for grid
  dd <- 11
  qq <- nrow(GAMA)
  minx <- min(RS[, 1])
  maxx <- max(RS[, 1])
  miny <- min(RS[, 2])
  maxy <- max(RS[, 2])

  # Create the indi vector
  indi <- NULL
  for (i in 1:nrow(GAMA)) {
    indi <- c(indi, rep(i, ncol(GAMA)))
  }

  # Create grid matrix
  se1 <- seq((minx + (maxx - minx) / (2 * dd)), (maxx - (maxx - minx) / (2 * dd)), (maxx - minx) / dd)
  se2 <- seq((miny + (maxy - miny) / (2 * dd)), (maxy - (maxy - miny) / (2 * dd)), (maxy - miny) / dd)
  MCELL <- array(NA, dim = c(length(se1), 2))

  MTCELL <- NULL
  for (i in 1:length(se1)) {
    for (j in 1:length(se2)) {
      MCELL[j, ] <- t(as.matrix(c(se1[i], se2[j])))
    }
    MTCELL <- rbind(MTCELL, MCELL)
  }

  S<-complex(real=sites[,1],imaginary=sites[,2])
  stemp<-H(nrow(sites))%*%S ## Removendo translação
  Sz<-t(H(nrow(sites)))%*%stemp ## configuração com translação removida
  Sz<-cbind(as.matrix(Re(Sz)),as.matrix(Im(Sz)))

  # Combine grid with observed sites
  DTotal <- rbind(MTCELL, sites)
  STotal <- complex(real = DTotal[, 1], imaginary = DTotal[, 2])
  stemp1 <- H(nrow(DTotal)) %*% STotal  # Remove translation
  SzT <- t(H(nrow(DTotal))) %*% stemp1  # Translated configuration
  SzT <- cbind(as.matrix(Re(SzT)), as.matrix(Im(SzT)))

  Sz1 <- SzT[(nrow(MTCELL) + 1):nrow(SzT), ]
  Sz2temp <- SzT[1:nrow(MTCELL), ]
  Sz2 <- Sz2temp[INDI, ]

  DD1<-Sz-Sz1
  DDD<-cbind(as.matrix(rep(DD1[1,1],nrow(Sz2))),as.matrix(rep(DD1[1,2],nrow(Sz2))) )
  Sz2c<-DDD+Sz2

  VecSNO <- as.matrix(as.vector(Sz2c))
  Szz1 <- rbind(Sz2c, Sz)
  distancias <- as.matrix(dist(Szz1))
  Rd1 <- exp(-bd * distancias)
  RA12 <- as.matrix(Rd1[1:nrow(Sz2c), (nrow(Sz2c) + 1):ncol(Rd1)])

  # Instructions for ind1 and ind2 (kept as originally provided)
  ind1 <- seq(1, 2 * ncol(response) * qq, 2 * ncol(response))
  ind2 <- seq(2 * ncol(response), 2 * ncol(response) * qq, 2 * ncol(response))

  # Initialize result matrix
  MatvecNO <- NULL
  vecm <- rbind(matrix(Sz[, 1], byrow = TRUE, ncol = 1), matrix(Sz[, 2], byrow = TRUE, ncol = 1))
  vecmu <- rbind(matrix(Sz2c[, 1], byrow = TRUE, ncol = 1), matrix(Sz2c[, 2], byrow = TRUE, ncol = 1))

  # Loop through MCMC samples
  for (l in 1:nrow(MTheta)) {

    DD <- NULL
    for (jj in 1:qq) {
      DD <- cbind(DD, as.matrix(MDef[l, ind1[jj]:ind2[jj]]))
    }

    VecDNO <- NULL
    mm <- rbind(Sz, Sz2c)
    distancias <- as.matrix(dist(mm))
    RD <- exp(-bd * distancias)
    Rd <- RD[1:ncol(response), 1:ncol(response)]
    Rdug <- RD[(ncol(response) + 1):nrow(mm), 1:ncol(response)]
    I2 <- diag(1, 2)

    for (jj in 1:qq) {
      rtem <- vecmu + (I2 %x% (Rdug %*% solve(Rd))) %*% (DD[, jj] - vecm)
      VecDNO <- cbind(VecDNO, rtem)
    }

    # Betaj matrix
    Betaj <- t(matrix(MTheta[l, ],nrow(FT) , nrow(response),byrow=T))
    vecyNO <- NULL

    # Interpolation loop for each time step
    for (t in 1:nrow(response)) {

      Dzl <- rbind(matrix(VecDNO[, indi[t]], ncol = 2), matrix(DD[, indi[t]], ncol = 2))
      SIGMAl <- Msigmak[l, indi[t]] * ((1 - Mkappa[l, indi[t]]) * gCorr(MPhi[l, indi[t]], Dzl) + Mkappa[l, indi[t]] * diag(1, nrow(Dzl)))
      SIGM <- SIGMAl[(nrow(Sz2) + 1):ncol(SIGMAl), (nrow(Sz2) + 1):ncol(SIGMAl)]
      SIGMNO <- SIGMAl[1:nrow(Sz2), 1:nrow(Sz2)]
      SIGEst <- t(as.matrix(SIGMAl[(nrow(Sz2) + 1):ncol(SIGMAl), 1:nrow(Sz2)]))

      mut <- t(FT) %*% as.matrix(Betaj[t, ])
      muNOt <- FTNO %*% as.matrix(Betaj[t, ])
      A <- muNOt + SIGEst %*% solve(SIGM) %*% (as.matrix(response[t, ]) - as.matrix(mut))
      B <- SIGMNO - SIGEst %*% solve(SIGM) %*% t(SIGEst)
      yNO <- t(as.matrix(mvrnorm(1, A, B)))

      vecyNO <- cbind(vecyNO, yNO)
    }

    MatvecNO <- rbind(MatvecNO, vecyNO)
  }
  return(MatvecNO)

}
