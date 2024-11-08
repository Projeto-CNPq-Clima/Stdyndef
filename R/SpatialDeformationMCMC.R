#' Fitting Spatiotemporal Model with Dynamic Deformation for Nonstationary Covariance Structures
#'
#' `SpatialDeformationMCMC` adjusts the parameters of a spatiotemporal model with dynamic deformation using MCMC sampling.
#'
#' @param response A matrix of observed values for the response variable, with dimensions `Txn`, where `T` is the length of the time series and `n` is the number of locations.
#' @param FT A covariate matrix with dimensions `pxn`, where `p` is the number of coefficients and `n` is the number of locations for time `t = 1`.
#' @param MatFFT A covariate matrix with dimensions `Tx(pxn)`. Each row contains the vectorized elements of the `FT` matrix for each time `t`.
#' @param GT An evolution matrix with dimensions `pxp`.
#' @param sites A matrix `nx2` containing the geographic coordinates of the locations.
#' @param GAMA A matrix with dimensions `qxr` such that `q * r = T`. Each row contains the positions of the set `AT = {1, ..., T}`, defining the subset `GAMA_i = AT[GAMA[i,]]` for `i = 1, ..., q`. The `GAMA` matrix must guarantee a partition of `AT`, meaning that `GAMA_i ∩ GAMA_j` is empty for all `i ≠ j`, and the union of all `GAMA_i` equals `AT`.
#' @param prior A list containing the hyperparameters for the prior distributions of the model parameters:
#'   * `theta_0 ~ N(m0, C0)`, with defaults `m0 = 0` and `C0 = diag(1000)`.
#'   * `Psi ~ Inverse-Wishart(S0, n0)`, with defaults `S0 = diag(0.01)` and `n0 = 1`.
#'   * `sigma_{k}^2 | pi_{k} ~ Inv-Gamma(axi + bxi, bepsilon / (1 - pi_k) + aepsilon / pi_k)`, with defaults `aepsilon = bepsilon = axi = bxi = 0.001`.
#'   * `phi ~ Inv-Gamma(aaphi, bbphi)`, with defaults `aaphi = bbphi = 0.001`.
#' @param iteration a value
#' @param burnin a value
#' @param jump a value
#'
#' @return A list containing the following adjusted parameters from the MCMC fitting procedure:
#'
#' \describe{
#'   \item{MDef}{A matrix of dimensions `((iteration-burnin)/jump) x (n x 2 x q)`. This matrix stores the deformation samples generated during the MCMC procedure. Each row contains the vectorized deformation values for the i-th iteration.}
#'   \item{MDef1T}{A binary vector indicating acceptance (1) or rejection (0) of the `d` parameter proposals in the Metropolis-Hastings step.}
#'   \item{MTheta}{A matrix of dimensions `((iteration-burnin)/jump) x (T x p)`. Each row contains a vectorized sample of the posterior distribution for `Theta` at the i-th iteration.}
#'   \item{MTheta0}{A matrix of dimensions `((iteration-burnin)/jump) x p`. Each row holds a vectorized sample of the posterior distribution for `Theta0` from the i-th iteration.}
#'   \item{Mlambda}{A matrix of dimensions `((iteration-burnin)/jump) x (n x 2 x q)`. This matrix stores the samples of `Z` values generated in the MCMC procedure, with each row representing the vectorized values of `Z` for each iteration.}
#'   \item{MPsi}{A matrix of dimensions `((iteration-burnin)/jump) x (p x p)`. This matrix contains the vectorized posterior samples of `Psi` for each MCMC iteration.}
#'   \item{Msigmak}{A matrix of dimensions `((iteration-burnin)/jump) x q`. This matrix contains the posterior samples of `sigmak` for each iteration, where each row represents values for the i-th iteration.}
#'   \item{Mkappa}{A matrix of dimensions `((iteration-burnin)/jump) x q`. This matrix contains the posterior samples of `kappa` values generated in the MCMC procedure, with each row representing the i-th iteration.}
#'   \item{MPhi}{A matrix of dimensions `((iteration-burnin)/jump) x q`. This matrix stores the posterior samples of `phi` values generated during the MCMC procedure. Each row corresponds to an iteration.}
#'   \item{MkappaT1}{A binary matrix of dimensions `((iteration-burnin)/jump) x q` indicating acceptance (1) or rejection (0) of the proposals for `kappa` parameters in the Metropolis-Hastings step.}
#'   \item{MPhiT1}{A binary matrix of dimensions `((iteration-burnin)/jump) x q` indicating acceptance (1) or rejection (0) of the proposals for `phi` parameters in the Metropolis-Hastings step.}
#'   \item{Mvarphi}{A vector of posterior samples for the `varphi` parameter obtained during the MCMC procedure.}
#' }
#'
#' @details This function performs MCMC sampling to fit a nonstationary spatiotemporal model with dynamic deformation, allowing for complex covariance structures in spatiotemporal data.
#'
#' @examples
#' # Example usage:
#' # Model fitting with simulated data
#' # result <- SpatialDeformationMCMC(Response, FT, MatFFT, GT, Sites, GAMA, prior)
#' @export





SpatialDeformationMCMC <- function(
    response, FT, MatFFT, GT, sites, GAMA,
    prior = list(
      m0 = as.matrix(rep(0, nrow(FT))),
      C0 = diag(1000, nrow(FT)),
      S0 = diag(0.01, ncol(GT)),
      n0 = 1,
      aepsilon = rep(0.01, nrow(GAMA)),
      axi = rep(0.01, nrow(GAMA)),
      bepsilon = rep(0.01, nrow(GAMA)),
      bxi = rep(0.01, nrow(GAMA)),
      aaphi = rep(0.01, nrow(GAMA)),
      bbphi = rep(0.01, nrow(GAMA)),
      sigmadd = 0.51,
      sigmal = 0.41
    ), iteration, burnin, jump) {
  # hiperparametros da lista prior

  m0 <- prior$m0
  mvarphi <- 1 / 18
  C0 <- prior$ C0
  S0 <- prior$S0
  n0 <- prior$n0
  n0est <- n0 + nrow(response)
  asigmad <- prior$asigmad
  bsigmad <- prior$bsigmad
  asigmala <- prior$asigmala
  bsigmala <- prior$bsigmala
  aepsilon <- prior$aepsilon
  axi <- prior$axi
  bepsilon <- prior$bepsilon
  bxi <- prior$bxi
  aaphi <- prior$aaphi
  bbphi <- prior$bbphi
  u1phi <- rep(500, nrow(GAMA))
  sigmadd <- prior$sigmadd
  sigmal <- prior$sigmal
  ###

  qq <- nrow(GAMA)
  Phis <- rep(0.51, qq)
  sigmas <- rep(0.51, qq)
  kappas <- rep(0.51, qq)
  varphi <- 0.51
  lambdas <- matrix(0.000, qq, ncol(response) * 2)

  # Valores computacionais
  uu1 <- rep(3.678794e-07, nrow(lambdas))
  lambdaant <- matrix(0, qq, ncol(lambdas))
  VV <- matrix(NA, qq, ncol(response)^2)
  RR <- matrix(NA, qq, ncol(response)^2)

  indi <- NULL

  for (i in 1:nrow(GAMA)) {
    indi <- c(indi, rep(i, ncol(GAMA)))
  }

  uuu1 <- rep(50, nrow(GAMA))
  uuu2 <- rep(50, nrow(GAMA))

  S <- complex(real = sites[, 1], imaginary = sites[, 2])
  stemp <- H(nrow(sites)) %*% S ## Removendo translação
  Sz <- t(H(nrow(sites))) %*% stemp ## configuração com translação removida
  Sz <- cbind(as.matrix(Re(Sz)), as.matrix(Im(Sz)))
  Z0 <- as.vector(Sz)

  MZ <- NULL
  for (l in 1:qq) {
    MZ <- rbind(MZ, t(as.matrix(as.vector(Sz))))
  }
  Zest <- matrix(NA, nrow(MZ), ncol(MZ))

  bd <- -2 * log(0.5) / (max(dist(Sz)))
  Delta <- gCorr(bd, Sz)

  Thetas <- matrix(0, nrow(response), ncol(GT))
  Theta0 <- as.matrix(rep(0, ncol(Thetas)))
  Psis <- diag(0.5, ncol(Thetas))

  MDef <- NULL
  MDefT <- NULL
  MDef1T <- NULL
  MTheta <- NULL
  MTheta0 <- NULL
  Mlambda <- NULL
  MPsi <- NULL
  Msigmak <- NULL
  Mkappa <- NULL
  MkappaT <- NULL
  MkappaT1 <- NULL
  MPhi <- NULL
  MPhiT <- NULL
  MPhiT1 <- NULL
  Mvarphi <- NULL

  for (j in 1:iteration) {
    if (j <= burnin) {
      DDefT <- NULL
      ttt <- NULL
      resphiT <- NULL


      for (s in 1:qq) {
        DEF <- matrix(MZ[s, ], ncol(response), 2)
        VV[s, ] <- as.vector(sigmas[s] * ((1 - kappas[s]) * gCorr(Phis[s], DEF) + kappas[s] * diag(1, ncol(response))))
        RR[s, ] <- as.vector(((1 - kappas[s]) * gCorr(Phis[s], DEF) + kappas[s] * diag(1, ncol(response))))
        Zest[s, ] <- MZ[s, ] - Z0

        temp <- amostrard(as.matrix(MZ[s, ]), Sz, Delta, MatFFT, response, Thetas, GAMA[s, ], kappas[s], sigmas[s], Phis[s], as.matrix(lambdas[s, ]), sigmadd, Delta, uu1[s])
        MZ[s, ] <- temp[[1]]
        DDefT <- c(DDefT, temp[[2]])

        if ((j %% 50) == 0) {
          uuu1[s] <- sintonizar(0.60, uuu1[s], c(MkappaT[, s], temp[[2]]), j)
          uuu2[s] <- uuu1[s]
          uu1[s] <- sintonizarN(0.20, uu1[s], c(MDefT[, s], temp[[2]]), j)
          u1phi[s] <- sintonizar(0.60, u1phi[s], c(MPhiT[, s], temp[[2]]), j)
        } else {

        }


        bbsigmak <- bxi[s] / (1 - kappas[s]) + bepsilon[s] / kappas[s] + 0.5 * SSUMmed(MatFFT[GAMA[s, ], ], response[GAMA[s, ], ], Thetas[GAMA[s, ], ], matrix(RR[s, ], ncol(response), ncol(response)))
        absigmak <- aepsilon[s] + axi + ncol(GAMA) * ncol(response) / 2
        sigmas[s] <- 1 / rgamma(1, absigmak, bbsigmak)


        temp <- tryCatch(amostrarkapa(MatFFT[GAMA[s, ], ], response[GAMA[s, ], ], Thetas[GAMA[s, ], ], kappas[s], sigmas[s], Phis[s], MZ[s, ], aepsilon[s], axi[s], bepsilon[s], bxi[s], uuu1[s], uuu2[s]), error = function(e) e)

        if (is.list(temp) == T) {
          temp <- amostrarkapa(MatFFT[GAMA[s, ], ], response[GAMA[s, ], ], Thetas[GAMA[s, ], ], kappas[s], sigmas[s], Phis[s], MZ[s, ], aepsilon[s], axi[s], bepsilon[s], bxi[s], uuu1[s], uuu2[s])
          kappas[s] <- temp[[1]]
          ttt <- c(ttt, temp[[2]])
        } else {
          ttt <- c(ttt, temp[[2]])
        }

        temp <- tryCatch(amostrarphi(MatFFT[GAMA[s, ], ], response[GAMA[s, ], ], Thetas[GAMA[s, ], ], kappas[s], sigmas[s], Phis[s], MZ[s, ], aaphi[s], bbphi[s], u1phi[s]), error = function(e) e)

        if (is.list(temp) == T) {
          temp <- amostrarphi(MatFFT[GAMA[s, ], ], response[GAMA[s, ], ], Thetas[GAMA[s, ], ], kappas[s], sigmas[s], Phis[s], MZ[s, ], aaphi[s], bbphi[s], u1phi[s])
          Phis[s] <- temp[[1]]
          resphiT <- c(resphiT, temp[[2]])
        } else {
          resphiT <- c(resphiT, 0)
        }
      }

      temp <- FFBSmod(MatFFT, GT, VV, Psis, m0, C0, response, Thetas, indi)
      Thetas <- temp

      CC0 <- solve(solve(C0) + t(GT) %*% solve(Psis) %*% GT)
      mm0 <- CC0 %*% (solve(C0) %*% m0 + t(GT) %*% solve(Psis) %*% as.matrix(Thetas[1, ]))

      Theta0 <- as.matrix(MASS::mvrnorm(1, mm0, CC0))

      MDefT <- rbind(MDefT, t(as.matrix(DDefT)))
      MkappaT <- rbind(MkappaT, t(as.matrix(ttt)))
      MPhiT <- rbind(MPhiT, t(as.matrix(resphiT)))


      temp <- FFBSdef(diag(1, 2 * ncol(response)), varphi * diag(1, 2 * ncol(response)), sigmadd * diag(1, 2) %x% Delta, sigmal * diag(1, 2) %x% Delta, as.matrix(rep(0, ncol(response) * 2)), diag(0.1, 2 * ncol(response)), Zest, lambdas)
      lambdas <- temp

      Thetaant <- rbind(t(Theta0), Thetas[1:(nrow(Thetas) - 1), ])
      S0est <- solve(S0 + SSUMPsi(GT, Thetas, Thetaant))
      temp <- solve(array(as.vector(rWishart(1, n0est, S0est)), dim = c(ncol(GT), ncol(GT))))
      Psis <- temp


      lambdaant[2:qq, ] <- lambdas[1:(qq - 1), ]


      SMA1 <- 0
      SMA2 <- 0
      for (kk in 1:qq) {
        rres1 <- t(as.matrix(lambdaant[kk, ])) %*% solve((sigmal * (diag(1, 2) %x% Delta))) %*% as.matrix(lambdas[kk, ])
        SMA1 <- SMA1 + rres1
        rres2 <- t(as.matrix(lambdaant[kk, ])) %*% solve((sigmal * (diag(1, 2) %x% Delta))) %*% as.matrix(lambdaant[kk, ])
        SMA2 <- SMA2 + rres2
      }

      VVarvarphi <- (SMA2 + mvarphi)^{
        -1
      }
      MMeanvarphi <- SMA1 * VVarvarphi
      varphi <- rnorm(1, MMeanvarphi, sqrt(VVarvarphi))
    } else {
      if (j == (burnin + 1)) {
        rm(MDefT, MkappaT, MPhiT)
      } else {

      }

      if ((j %% jump) == 0) {
        rrr <- NULL
        DDef <- NULL
        DDef1T <- NULL
        rrr1 <- NULL
        resphi <- NULL
        reskappa <- NULL
        resphi <- NULL
        resphiT <- NULL

        for (s in 1:qq) {
          DEF <- matrix(MZ[s, ], ncol(response), 2)
          VV[s, ] <- as.vector(sigmas[s] * ((1 - kappas[s]) * gCorr(Phis[s], DEF) + kappas[s] * diag(1, ncol(response))))
          RR[s, ] <- as.vector(((1 - kappas[s]) * gCorr(Phis[s], DEF) + kappas[s] * diag(1, ncol(response))))
          Zest[s, ] <- MZ[s, ] - Z0

          temp <- amostrard(as.matrix(MZ[s, ]), Sz, Delta, MatFFT, response, Thetas, GAMA[s, ], kappas[s], sigmas[s], Phis[s], as.matrix(lambdas[s, ]), sigmadd, Delta, uu1[s])
          MZ[s, ] <- temp[[1]]
          DDef <- cbind(DDef, t(temp[[1]]))
          DDef1T <- c(DDef1T, temp[[2]])

          bbsigmak <- bxi[s] / (1 - kappas[s]) + bepsilon[s] / kappas[s] + 0.5 * SSUMmed(MatFFT[GAMA[s, ], ], response[GAMA[s, ], ], Thetas[GAMA[s, ], ], matrix(RR[s, ], ncol(response), ncol(response)))
          absigmak <- aepsilon[s] + axi + ncol(GAMA) * ncol(response) / 2
          sigmas[s] <- 1 / rgamma(1, absigmak, bbsigmak)
          rrr <- c(rrr, sigmas[s])

          temp <- tryCatch(amostrarkapa(MatFFT[GAMA[s, ], ], response[GAMA[s, ], ], Thetas[GAMA[s, ], ], kappas[s], sigmas[s], Phis[s], MZ[s, ], aepsilon[s], axi[s], bepsilon[s], bxi[s], uuu1[s], uuu2[s]), error = function(e) e)

          if (is.list(temp) == T) {
            temp <- amostrarkapa(MatFFT[GAMA[s, ], ], response[GAMA[s, ], ], Thetas[GAMA[s, ], ], kappas[s], sigmas[s], Phis[s], MZ[s, ], aepsilon[s], axi[s], bepsilon[s], bxi[s], uuu1[s], uuu2[s])
            kappas[s] <- temp[[1]]
            rrr1 <- c(rrr1, kappas[s])
            reskappa <- c(reskappa, temp[[2]])
          } else {
            rrr1 <- c(rrr1, kappas[s])
            reskappa <- c(reskappa, 0)
          }

          temp <- tryCatch(amostrarphi(MatFFT[GAMA[s, ], ], response[GAMA[s, ], ], Thetas[GAMA[s, ], ], kappas[s], sigmas[s], Phis[s], MZ[s, ], aaphi[s], bbphi[s], u1phi[s]), error = function(e) e)

          if (is.list(temp) == T) {
            temp <- amostrarphi(MatFFT[GAMA[s, ], ], response[GAMA[s, ], ], Thetas[GAMA[s, ], ], kappas[s], sigmas[s], Phis[s], MZ[s, ], aaphi[s], bbphi[s], u1phi[s])
            Phis[s] <- temp[[1]]
            resphiT <- c(resphiT, temp[[2]])
            resphi <- c(resphi, temp[[1]])
          } else {
            resphiT <- c(resphiT, 0)
          }
        }

        Msigmak <- rbind(Msigmak, t(as.matrix(rrr)))

        Mkappa <- rbind(Mkappa, t(as.matrix(rrr1)))
        MkappaT1 <- rbind(MkappaT1, t(as.matrix(reskappa)))

        MPhi <- rbind(MPhi, t(as.matrix(resphi)))
        MPhiT1 <- rbind(MPhiT1, t(as.matrix(resphiT)))

        MDef <- rbind(MDef, DDef)
        MDef1T <- rbind(MDef1T, t(as.matrix(DDef1T)))

        temp <- FFBSmod(MatFFT, GT, VV, Psis, m0, C0, response, Thetas, indi)
        Thetas <- temp
        MTheta <- rbind(MTheta, t(as.matrix(as.vector(Thetas))))

        CC0 <- solve(solve(C0) + t(GT) %*% solve(Psis) %*% GT)
        mm0 <- CC0 %*% (solve(C0) %*% m0 + t(GT) %*% solve(Psis) %*% as.matrix(Thetas[1, ]))

        Theta0 <- as.matrix(mvrnorm(1, mm0, CC0))
        MTheta0 <- rbind(MTheta0, t(Theta0))

        temp <- FFBSdef(diag(1, 2 * ncol(response)), varphi * diag(1, 2 * ncol(response)), sigmadd * diag(1, 2) %x% Delta, sigmal * diag(1, 2) %x% Delta, as.matrix(rep(0, ncol(response) * 2)), diag(0.1, 2 * ncol(response)), Zest, lambdas)
        lambdas <- temp
        Mlambda <- rbind(Mlambda, t(as.matrix(as.vector(lambdas))))

        Thetaant <- rbind(t(Theta0), Thetas[1:(nrow(Thetas) - 1), ])
        S0est <- solve(S0 + SSUMPsi(GT, Thetas, Thetaant))
        temp <- solve(array(as.vector(rWishart(1, n0est, S0est)), dim = c(ncol(GT), ncol(GT))))
        Psis <- temp
        MPsi <- rbind(MPsi, t(as.matrix(as.vector(Psis))))

        lambdaant[2:qq, ] <- lambdas[1:(qq - 1), ]

        SMA1 <- 0
        SMA2 <- 0
        for (kk in 1:qq) {
          rres1 <- t(as.matrix(lambdaant[kk, ])) %*% solve((sigmal * (diag(1, 2) %x% Delta))) %*% as.matrix(lambdas[kk, ])
          SMA1 <- SMA1 + rres1
          rres2 <- t(as.matrix(lambdaant[kk, ])) %*% solve((sigmal * (diag(1, 2) %x% Delta))) %*% as.matrix(lambdaant[kk, ])
          SMA2 <- SMA2 + rres2
        }

        VVarvarphi <- (SMA2 + mvarphi)^{
          -1
        }
        MMeanvarphi <- SMA1 * VVarvarphi
        varphi <- rnorm(1, MMeanvarphi, sqrt(VVarvarphi))
        Mvarphi <- c(Mvarphi, varphi)
      } else {
        rrr <- NULL
        DDef <- NULL
        rrr1 <- NULL
        resphi <- NULL

        for (s in 1:qq) {
          DEF <- matrix(MZ[s, ], ncol(response), 2)
          VV[s, ] <- as.vector(sigmas[s] * ((1 - kappas[s]) * gCorr(Phis[s], DEF) + kappas[s] * diag(1, ncol(response))))
          RR[s, ] <- as.vector(((1 - kappas[s]) * gCorr(Phis[s], DEF) + kappas[s] * diag(1, ncol(response))))
          Zest[s, ] <- MZ[s, ] - Z0

          temp <- amostrard(as.matrix(MZ[s, ]), Sz, Delta, MatFFT, response, Thetas, GAMA[s, ], kappas[s], sigmas[s], Phis[s], as.matrix(lambdas[s, ]), sigmadd, Delta, uu1[s])
          MZ[s, ] <- temp[[1]]
          DDef <- cbind(DDef, t(temp[[1]]))

          bbsigmak <- bxi[s] / (1 - kappas[s]) + bepsilon[s] / kappas[s] + 0.5 * SSUMmed(MatFFT[GAMA[s, ], ], response[GAMA[s, ], ], Thetas[GAMA[s, ], ], matrix(RR[s, ], ncol(response), ncol(response)))
          absigmak <- aepsilon[s] + axi + ncol(GAMA) * ncol(response) / 2
          sigmas[s] <- 1 / rgamma(1, absigmak, bbsigmak)
          rrr <- c(rrr, sigmas[s])

          temp <- tryCatch(amostrarkapa(MatFFT[GAMA[s, ], ], response[GAMA[s, ], ], Thetas[GAMA[s, ], ], kappas[s], sigmas[s], Phis[s], MZ[s, ], aepsilon[s], axi[s], bepsilon[s], bxi[s], uuu1[s], uuu2[s]), error = function(e) e)

          if (is.list(temp) == T) {
            temp <- amostrarkapa(MatFFT[GAMA[s, ], ], response[GAMA[s, ], ], Thetas[GAMA[s, ], ], kappas[s], sigmas[s], Phis[s], MZ[s, ], aepsilon[s], axi[s], bepsilon[s], bxi[s], uuu1[s], uuu2[s])
            kappas[s] <- temp[[1]]
            rrr1 <- c(rrr1, kappas[s])
          } else {
            rrr1 <- c(rrr1, kappas[s])
          }

          temp <- tryCatch(amostrarphi(MatFFT[GAMA[s, ], ], response[GAMA[s, ], ], Thetas[GAMA[s, ], ], kappas[s], sigmas[s], Phis[s], MZ[s, ], aaphi[s], bbphi[s], u1phi[s]), error = function(e) e)

          if (is.list(temp) == T) {
            temp <- amostrarphi(MatFFT[GAMA[s, ], ], response[GAMA[s, ], ], Thetas[GAMA[s, ], ], kappas[s], sigmas[s], Phis[s], MZ[s, ], aaphi[s], bbphi[s], u1phi[s])
            Phis[s] <- temp[[1]]
            resphi <- c(resphi, temp[[1]])
          } else {
            resphi <- c(resphi, temp[[1]])
          }
        }

        temp <- FFBSmod(MatFFT, GT, VV, Psis, m0, C0, response, Thetas, indi)
        Thetas <- temp

        CC0 <- solve(solve(C0) + t(GT) %*% solve(Psis) %*% GT)
        mm0 <- CC0 %*% (solve(C0) %*% m0 + t(GT) %*% solve(Psis) %*% as.matrix(Thetas[1, ]))

        Theta0 <- as.matrix(mvrnorm(1, mm0, CC0))

        temp <- FFBSdef(diag(1, 2 * ncol(response)), varphi * diag(1, 2 * ncol(response)), sigmadd * diag(1, 2) %x% Delta, sigmal * diag(1, 2) %x% Delta, as.matrix(rep(0, ncol(response) * 2)), diag(0.1, 2 * ncol(response)), Zest, lambdas)
        lambdas <- temp

        Thetaant <- rbind(t(Theta0), Thetas[1:(nrow(Thetas) - 1), ])
        S0est <- solve(S0 + SSUMPsi(GT, Thetas, Thetaant))
        temp <- solve(array(as.vector(rWishart(1, n0est, S0est)), dim = c(ncol(GT), ncol(GT))))
        Psis <- temp


        lambdaant[2:qq, ] <- lambdas[1:(qq - 1), ]


        SMA1 <- 0
        SMA2 <- 0
        for (kk in 1:qq) {
          rres1 <- t(as.matrix(lambdaant[kk, ])) %*% solve((sigmal * (diag(1, 2) %x% Delta))) %*% as.matrix(lambdas[kk, ])
          SMA1 <- SMA1 + rres1
          rres2 <- t(as.matrix(lambdaant[kk, ])) %*% solve((sigmal * (diag(1, 2) %x% Delta))) %*% as.matrix(lambdaant[kk, ])
          SMA2 <- SMA2 + rres2
        }

        VVarvarphi <- (SMA2 + mvarphi)^{
          -1
        }
        MMeanvarphi <- SMA1 * VVarvarphi
        varphi <- rnorm(1, MMeanvarphi, sqrt(VVarvarphi))
      }
    }

    print(j)
  }

  resultss <- list(MDef,
                               MDef1T,
                               MTheta,
                               MTheta0,
                               Mlambda,
                               MPsi,
                               Msigmak,
                               Mkappa,
                               MkappaT1,
                               MPhi,
                               MPhiT1,
                               Mvarphi,
                               bd)


  names(resultss) <- c(
    "MDef",
    "MDef1T",
    "MTheta",
    "MTheta0",
    "Mlambda",
    "MPsi",
    "Msigmak",
    "Mkappa",
    "MkappaT1",
    "MPhi",
    "MPhiT1",
    "Mvarphi",
    "bd"
  )


  return(resultss)
}
