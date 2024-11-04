sintonizarN <- function(taxa, tau, mat, i) {
  mat <- as.matrix(mat)


  mater <- (1 / 50) * sum(mat[(i - 49):i, 1])

  if (mater >= taxa) {
    delta <- min(0.01, (i / 50 + 1)^(-0.5))
    temp4 <- log(tau) + delta
    temp5 <- exp(temp4)
    return(temp5)
  } else {
    delta <- min(0.01, (i / 50 + 1)^(-0.5))
    temp4 <- log(tau) - delta
    temp5 <- exp(temp4)
    return(temp5)
  }
}

#

sintonizar <- function(taxa, tau, mat, i) {
  mat <- as.matrix(mat)



  mater <- (1 / 50) * sum(mat[(i - 49):i, 1])

  if (mater >= taxa) {
    delta <- min(0.01, (i / 50 + 1)^(-0.5))
    temp4 <- log(tau) - delta
    temp5 <- exp(temp4)
    return(temp5)
  } else {
    delta <- min(0.01, (i / 50 + 1)^(-0.5))
    temp4 <- log(tau) + delta
    temp5 <- exp(temp4)
    return(temp5)
  }
}
#
FFBSmod <- function(MFFt, GGt, MVVt, WWt, mm0, CC0, yy, ttheta, indices) {
  tempA <- matrix(NA, nrow(ttheta), ncol(ttheta))
  Tt <- nrow(yy)
  mm <- matrix(0, nrow(yy), nrow(mm0))
  CC <- matrix(0, nrow(yy), (nrow(CC0) * ncol(CC0)))

  mmant <- mm0
  CCant <- CC0

  for (i in 1:nrow(yy)) {
    aat <- GGt %*% mmant
    RRt <- GGt %*% CCant %*% t(GGt) + WWt
    FFt <- matrix(MFFt[i, ], ncol(ttheta), ncol(yy))
    fft <- t(FFt) %*% aat
    VVt <- matrix(MVVt[indices[i], ], ncol(yy), ncol(yy))
    QQt <- t(FFt) %*% RRt %*% FFt + VVt
    AAt <- RRt %*% FFt %*% solve(QQt)
    eet <- as.matrix(yy[i, ]) - fft
    mmt <- aat + AAt %*% eet
    mm[i, ] <- mmt
    CCt <- RRt - AAt %*% QQt %*% t(AAt)
    CC[i, ] <- as.vector(CCt)
    mmant <- mmt
    CCant <- CCt
  }

  TTt <- nrow(yy)
  mmm <- as.matrix(mm[TTt, ])
  CCC <- matrix(CC[TTt, ], ncol(mm), ncol(mm))
  tempA[TTt, ] <- as.matrix(mvrnorm(1, mmm, CCC))

  for (j in (TTt - 1):1) {
    mmm <- as.matrix(mm[j, ])
    CCC <- matrix(CC[j, ], ncol(mm), ncol(mm))
    varA <- solve(t(GGt) %*% solve(WWt) %*% GGt + solve(CCC))
    mediaA <- varA %*% (t(GGt) %*% solve(WWt) %*% as.matrix(tempA[(j + 1), ]) + solve(CCC) %*% mmm)
    tempA[j, ] <- as.matrix(mvrnorm(1, mediaA, varA))
  }
  tempA
}
#
FFBSdef <- function(FFt, GGt, VVt, WWt, mm0, CC0, yy, ttheta) {
  tempA <- matrix(NA, nrow(ttheta), ncol(ttheta))
  Tt <- nrow(yy)
  mm <- matrix(0, nrow(yy), nrow(mm0))
  CC <- matrix(0, nrow(yy), (nrow(CC0) * ncol(CC0)))

  mmant <- mm0
  CCant <- CC0

  for (i in 1:nrow(yy)) {
    aat <- GGt %*% mmant
    RRt <- GGt %*% CCant %*% t(GGt) + WWt
    fft <- t(FFt) %*% aat
    QQt <- t(FFt) %*% RRt %*% FFt + VVt
    AAt <- RRt %*% FFt %*% solve(QQt)
    eet <- as.matrix(yy[i, ]) - fft
    mmt <- aat + AAt %*% eet
    mm[i, ] <- mmt
    CCt <- RRt - AAt %*% QQt %*% t(AAt)
    CC[i, ] <- as.vector(CCt)
    mmant <- mmt
    CCant <- CCt
  }

  TTt <- nrow(yy)
  mmm <- as.matrix(mm[TTt, ])
  CCC <- matrix(CC[TTt, ], ncol(mm), ncol(mm))
  tempA[TTt, ] <- as.matrix(mvrnorm(1, mmm, CCC))

  for (j in (TTt - 1):1) {
    mmm <- as.matrix(mm[j, ])
    CCC <- matrix(CC[j, ], ncol(mm), ncol(mm))
    varA <- solve(t(GGt) %*% solve(WWt) %*% GGt + solve(CCC))
    mediaA <- varA %*% (t(GGt) %*% solve(WWt) %*% as.matrix(tempA[(j + 1), ]) + solve(CCC) %*% mmm)
    tempA[j, ] <- as.matrix(mvrnorm(1, mediaA, varA))
  }

  tempA
}
#
SSUMPsi <- function(MAT, THETA, THETAant) {
  suma <- matrix(0, ncol(THETA), ncol(THETA))

  for (j in 1:nrow(THETA)) {
    res <- (as.matrix(THETA[j, ]) - MAT %*% as.matrix(THETAant[j, ])) %*% t((as.matrix(THETA[j, ]) - MAT %*% as.matrix(THETAant[j, ])))
    suma <- suma + res
  }

  suma
}
#
SSUMd <- function(MAT, THETA, THETAant, DDD) {
  suma <- 0

  for (j in 1:nrow(THETA)) {
    res <- t(as.matrix(THETA[j, ]) - MAT %*% as.matrix(THETAant[j, ])) %*% solve(DDD) %*% (as.matrix(THETA[j, ]) - MAT %*% as.matrix(THETAant[j, ]))
    suma <- suma + res
  }

  suma
}
#
SSUMmed <- function(MAT, YY, THETA, DDD) {
  suma <- 0

  for (k in 1:nrow(YY)) {
    res <- t((as.matrix(YY[k, ]) - t(matrix(MAT[k, ], ncol(THETA), ncol(YY))) %*% as.matrix(THETA[k, ]))) %*% solve(DDD) %*% (as.matrix(YY[k, ]) - t(matrix(MAT[k, ], ncol(THETA), ncol(YY))) %*% as.matrix(THETA[k, ]))
    suma <- suma + res
  }

  suma
}
#
confh <- function(vetor) {
  n <- length(vetor)

  w <- as.matrix(H(n) %*% vetor)
  escala <- as.numeric(sqrt(t(Conj(w)) %*% w))
  z <- (1 / escala) * w
  res <- list(z, escala)
}
#
H <- function(k) {
  matriz <- array(0, dim = c((k - 1), k))

  for (i in 1:(k - 1)) {
    for (j in 1:k) {
      if (j <= i) {
        matriz[i, j] <- -(i * (i + 1))^(-1 / 2)
      } else {
        if (j == (i + 1)) {
          matriz[i, j] <- -i * (-(i * (i + 1))^(-1 / 2))
        } else {

        }
      }
    }
  }

  matriz
}
#
gCorr <- function(bb, dd) {
  n <- nrow(dd)
  R <- exp(-bb * (as.matrix(dist(dd))))
  R
}
#
amostrard <- function(DDef, Sr, Rde, MFFT, yTT, TTheta, indGama, kkappa, sssigma, PPhi, llambda, ssigmad, DDelta, u1) {
  QQ <- ssigmad * diag(1, 2) %x% DDelta
  Z0 <- as.matrix(as.vector(Sr))
  n <- nrow(DDef) / 2
  Do <- complex(real = DDef[1:(n), 1], imaginary = DDef[(n + 1):(2 * n), 1]) # Coordenadas complexas
  ztemp <- confh(Do) ## Removendo translação e escala
  z <- ztemp[[1]] ## Pre-forma de D actual
  esc <- ztemp[[2]] ## valor da escala removida
  marcoc1 <- complex(real = Sr[, 1], imaginary = Sr[, 2])
  Zotemp <- confh(marcoc1)
  Zo <- Zotemp[[1]] # Pre-forma de referencia
  theta <- Arg(t(Conj(z)) %*% Zo)
  ww <- as.vector(exp(theta * (0 + 1i))) * z
  zd <- t(H(n)) %*% ww ## Configuração Pre-forma de D actual
  vectD <- c(Re(zd), Im(zd)) # Configuração pre-forma de D actual vetorizada
  ssigma <- u1 * diag(rep(1, 2)) %x% Rde # Matriz de covariancia proposta função de transição
  vectdprop <- mvrnorm(1, vectD, ssigma) # Geração de configuração Pre-forma proposta vetorizada
  Dp <- cbind(as.matrix(vectdprop[1:n]), as.matrix(vectdprop[(n + 1):(2 * n)])) # Configuração proposta pre-forma
  dprop <- esc * Dp # D proposto escala original
  vecdprop <- as.matrix(c(dprop[, 1], dprop[, 2]))

  INN <- diag(1, ncol(yTT))

  ddeff <- cbind(as.matrix(DDef[1:(nrow(DDef) / 2), 1]), as.matrix(DDef[(nrow(DDef) / 2 + 1):nrow(DDef), 1]))
  SSS <- sssigma * ((1 - kkappa) * gCorr(PPhi, ddeff) + kkappa * INN)

  SUM <- 0
  for (j in indGama) {
    FFT <- matrix(MFFT[j, ], ncol(TTheta), ncol(yTT))
    SUM <- SUM + t(as.matrix(yTT[j, ]) - t(FFT) %*% as.matrix(TTheta[j, ])) %*% solve(SSS) %*% (as.matrix(yTT[j, ]) - t(FFT) %*% as.matrix(TTheta[j, ]))
  }

  logp <- -(length(indGama) / 2) * log(det(SSS)) - 0.5 * SUM - 0.5 * t(DDef - Z0 - llambda) %*% solve(QQ) %*% (DDef - Z0 - llambda)

  ddeffprop <- cbind(as.matrix(vecdprop[1:(nrow(vecdprop) / 2), 1]), as.matrix(vecdprop[(nrow(vecdprop) / 2 + 1):nrow(vecdprop), 1]))
  SSSp <- sssigma * ((1 - kkappa) * gCorr(PPhi, ddeffprop) + kkappa * INN)

  SUMp <- 0
  for (j in indGama) {
    FFT <- matrix(MFFT[j, ], ncol(TTheta), ncol(yTT))
    SUMp <- SUMp + t(as.matrix(yTT[j, ]) - t(FFT) %*% as.matrix(TTheta[j, ])) %*% solve(SSSp) %*% (as.matrix(yTT[j, ]) - t(FFT) %*% as.matrix(TTheta[j, ]))
  }

  logpprop <- -(length(indGama) / 2) * log(det(SSSp)) - 0.5 * SUMp - 0.5 * t(vecdprop - Z0 - llambda) %*% solve(QQ) %*% (vecdprop - Z0 - llambda)

  prob <- min(exp((logpprop) - (logp)), 1)

  u <- runif(1, 0, 1)

  if (u < prob) {
    dprox <- vecdprop

    rejei <- 1
  } else {
    dprox <- DDef
    rejei <- 0
  }






  res <- list(dprox, rejei)
  res
}
#
amostrarkapa <- function(MFFT, yTT, TTheta, kkappa, ssigma, PPhi, MDEF, aa1z, aa1e, bb1z, bb1e, u1, u2) {
  INN <- diag(1, ncol(yTT))
  ddeff <- matrix(MDEF, ncol(yTT), 2)
  SSS <- ssigma * ((1 - kkappa) * gCorr(PPhi, ddeff) + kkappa * INN)
  SUM <- SSUMmed(MFFT, yTT, TTheta, SSS)

  logp <- (aa1z - 1) * log(kkappa) + (aa1e - 1) * log(1 - kkappa) - (aa1e + aa1z) * log(bb1z * kkappa + bb1e * (1 - kkappa)) - (nrow(MFFT) / 2) * log(det(SSS)) - 0.5 * SUM

  kkappaprop <- rbeta(1, kkappa * (u1 + u2), (u1 + u2) * (1 - kkappa))

  if ((kkappaprop <= 0.0001) | (kkappaprop >= 0.9999)) {
    return(list(kkappa, 0))
  } else {

  }

  SSSp <- ssigma * ((1 - kkappaprop) * gCorr(PPhi, ddeff) + kkappaprop * INN)
  SUMp <- SSUMmed(MFFT, yTT, TTheta, SSSp)

  logpprop <- (aa1z - 1) * log(kkappaprop) + (aa1e - 1) * log(1 - kkappaprop) - (aa1e + aa1z) * log(bb1z * kkappaprop + bb1e * (1 - kkappaprop)) - (nrow(MFFT) / 2) * log(det(SSSp)) - 0.5 * SUMp

  logprob <- logpprop + log(dbeta(kkappa, kkappaprop * (u1 + u2), (u1 + u2) * (1 - kkappaprop))) - (logp + log(dbeta(kkappaprop, kkappa * (u1 + u2), (u1 + u2) * (1 - kkappa))))

  prob <- min(c(1, exp(logprob)))

  u <- runif(1, 0, 1)


  prob <- min(c(1, exp(logprob)))

  u <- runif(1, 0, 1)

  if (u < prob) {
    kapaprox <- kkappaprop

    rejei <- 1
  } else {
    kapaprox <- kkappa

    rejei <- 0
  }



  res <- list(kapaprox, rejei)
  res
}
#
amostrarphi <- function(MFFT, yTT, TTheta, kkappa, ssigma, PPhi, MDEF, aa1, bb1, u1) {
  INN <- diag(1, ncol(yTT))
  ddeff <- matrix(MDEF, ncol(yTT), 2)
  SSS <- ssigma * ((1 - kkappa) * gCorr(PPhi, ddeff) + kkappa * INN)
  SUM <- SSUMmed(MFFT, yTT, TTheta, SSS)

  logp <- -(nrow(MFFT) / 2) * log(det(SSS)) - 0.5 * SUM - bb1 / PPhi - (aa1 + 1) * log(PPhi)

  PPhiprop <- rgamma(1, shape = PPhi * u1, rate = u1)
  SSSp <- ssigma * ((1 - kkappa) * gCorr(PPhiprop, ddeff) + kkappa * INN)
  SUMp <- SSUMmed(MFFT, yTT, TTheta, SSSp)

  logpprop <- -(nrow(MFFT) / 2) * log(det(SSSp)) - 0.5 * SUMp - bb1 / PPhiprop - (aa1 + 1) * log(PPhiprop)

  logprob <- logpprop + log(dgamma(PPhi, shape = PPhiprop * u1, rate = u1)) - (logp + log(dgamma(PPhiprop, shape = PPhi * u1, rate = u1)))
  prob <- min(c(1, exp(logprob)))


  u <- runif(1, 0, 1)

  if (u < prob) {
    bprox <- PPhiprop

    rejei <- 1
  } else {
    bprox <- PPhi

    rejei <- 0
  }



  res <- list(bprox, rejei)
  res
}



# Define the logNor function
logNor <- function(MFFT, yTT, TTheta, kkappa, ssigma, PPhi, MDEF) {
  INN <- diag(1, ncol(yTT))
  ddeff <- matrix(MDEF, ncol(yTT), 2)
  SSS <- ssigma * ((1 - kkappa) * gCorr(PPhi, ddeff) + kkappa * INN)

  matres <- NULL
  for (j in 1:nrow(MFFT)) {
    mmm <- t(matrix(MFFT[j, ], ncol(TTheta), ncol(yTT))) %*% as.matrix(TTheta[j, ])
    matres <- c(matres, dmvnorm(yTT[j, ], mean = mmm, sigma = SSS, log = TRUE))
  }

  sum(matres)
}

# Define the logvero function
logvero <- function(MatT, Response, TTheta, kkappa, ssigma, PPhi, MDEF, IND) {
  tempsoma <- 0

  for (j in 1:nrow(IND)) {
    tempsoma <- tempsoma + logNor(MatT[IND[j, ], ], Response[IND[j, ], ], TTheta[IND[j, ], ],
                                  kkappa[j], ssigma[j], PPhi[j], MDEF[j, ])
  }

  tempsoma
}


