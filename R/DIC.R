#' Title
#'
#' @param output a value
#' @param response a value
#' @param FT a value
#' @param MatFFT a value
#' @param GT a value
#' @param sites a value
#' @param GAMA a value
#'
#' @return
#' @export
#'
#' @examples
#'
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
