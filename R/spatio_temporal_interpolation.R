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
#' @examples
#'library(graphics)
#'library(animation)
#'library(MASS)
#'library(Stdyndef)
#'
#'data(temperature)
#'data(FT)
#'data(MatFFT)
#'data(GT)
#'data(sites)
#'data(GAMA)
#'data(Map)
#'data(FTNO)
#'
#'Mod=SpatialDeformationMCMC(response=temperature,FT=FT,MatFFT=MatFFT,GT=GT,sites=sites,GAMA=GAMA,iteration=100,burnin=50,jump=1)
#'
#'INDI=c(4,15,16,26:28,36:39,46:66,69:77,80:88,92:99,105:108,119)
#'
#'Fig=(spatio_temporal_interpolation(Mod,temperature, FT, MatFFT, GT, sites, INDI, GAMA, FTNO, Map))
#'
#'
#'vecyest=apply(Fig,2,mean)
#'vecyp025=apply(Fig,2,quantile,prob=c(0.025))
#'vecyp0975=apply(Fig,2,quantile,prob=c(0.975))
#'
#'
#'yest=t(matrix(vecyest,length(INDI),nrow(temperature)))
#'YEAR=c("2007/Jan./Summer","2007/Feb./Summer","2007/Mar./Summer","2007/Apr./Autumn","2007/May/Autumn","2007/June/Autumn","2007/July/Winter","2007/Aug./Winter","2007/Sept./Winter","2007/Oct./Spring","2007/Nov./Spring","2007/Dec./Spring","2008/Jan./Summer","2008/Feb./Summer","2008/Mar./Summer","2008/Apr./Autumn","2008/May/Autumn","2008/June/Autumn","2008/July/Winter","2008/Aug./Winter","2008/Sept./Winter","2008/Oct./Spring","2008/Nov./Spring","2008/Dec./Spring",
#'       "2009/Jan./Summer","2009/Feb./Summer","2009/Mar./Summer","2009/Apr./Autumn","2009/May/Autumn","2009/June/Autumn","2009/July/Winter","2009/Aug./Winter","2009/Sept./Winter","2009/Oct./Spring","2009/Nov./Spring","2009/Dec./Spring",
#'       "2010/Jan./Summer","2010/Feb./Summer","2010/Mar./Summer","2010/Apr./Autumn","2010/May/Autumn","2010/June/Autumn","2010/July/Winter","2010/Aug./Winter","2010/Sept./Winter","2010/Oct./Spring","2010/Nov./Spring","2010/Dec./Spring",
#'       "2011/Jan./Summer","2011/Feb./Summer","2011/Mar./Summer","2011/Apr./Autumn","2011/May/Autumn","2011/June/Autumn","2011/July/Winter","2011/Aug./Winter","2011/Sept./Winter","2011/Oct./Spring","2011/Nov./Spring","2011/Dec./Spring",
#'       "2012/Jan./Summer","2012/Feb./Summer","2012/Mar./Summer","2012/Apr./Autumn","2012/May/Autumn","2012/June/Autumn","2012/July/Winter","2012/Aug./Winter","2012/Sept./Winter","2012/Oct./Spring","2012/Nov./Spring","2012/Dec./Spring",
#'       "2013/Jan./Summer","2013/Feb./Summer","2013/Mar./Summer","2013/Apr./Autumn","2013/May/Autumn","2013/June/Autumn","2013/July/Winter","2013/Aug./Winter","2013/Sept./Winter","2013/Oct./Spring","2013/Nov./Spring","2013/Dec./Spring",
#'       "2014/Jan./Summer","2014/Feb./Summer","2014/Mar./Summer","2014/Apr./Autumn","2014/May/Autumn","2014/June/Autumn","2014/July/Winter","2014/Aug./Winter","2014/Sept./Winter","2014/Oct./Spring","2014/Nov./Spring","2014/Dec./Spring",
#'       "2015/Jan./Summer","2015/Feb./Summer","2015/Mar./Summer","2015/Apr./Autumn","2015/May/Autumn","2015/June/Autumn","2015/July/Winter","2015/Aug./Winter","2015/Sept./Winter","2015/Oct./Spring","2015/Nov./Spring","2015/Dec./Spring",
#'       "2016/Jan./Summer","2016/Feb./Summer","2016/Mar./Summer","2016/Apr./Autumn","2016/May/Autumn","2016/June/Autumn","2016/July/Winter","2016/Aug./Winter","2016/Sept./Winter","2016/Oct./Spring","2016/Nov./Spring","2016/Dec./Spring",
#'       "2017/Jan./Summer","2017/Feb./Summer","2017/Mar./Summer","2017/Apr./Autumn","2017/May/Autumn","2017/June/Autumn","2017/July/Winter","2017/Aug./Winter","2017/Sept./Winter","2017/Oct./Spring","2017/Nov./Spring","2017/Dec./Spring",
#'       "2018/Jan./Summer","2018/Feb./Summer","2018/Mar./Summer","2018/Apr./Autumn","2018/May/Autumn","2018/June/Autumn","2018/July/Winter","2018/Aug./Winter","2018/Sept./Winter","2018/Oct./Spring","2018/Nov./Spring","2018/Dec./Spring",
#'       "2019/Jan./Summer","2019/Feb./Summer","2019/Mar./Summer","2019/Apr./Autumn","2019/May/Autumn","2019/June/Autumn","2019/July/Winter","2019/Aug./Winter","2019/Sept./Winter","2019/Oct./Spring","2019/Nov./Spring","2019/Dec./Spring",
#'       "2020/Jan./Summer","2020/Feb./Summer","2020/Mar./Summer","2020/Apr./Autumn","2020/May/Autumn","2020/June/Autumn","2020/July/Winter","2020/Aug./Winter","2020/Sept./Winter","2020/Oct./Spring","2020/Nov./Spring","2020/Dec./Spring",
#'       "2021/Jan./Summer","2021/Feb./Summer","2021/Mar./Summer","2021/Apr./Autumn","2021/May/Autumn","2021/June/Autumn","2021/July/Winter","2021/Aug./Winter","2021/Sept./Winter","2021/Oct./Spring","2021/Nov./Spring","2021/Dec./Spring")
#'infz=min(yest)#8.38#min(yest)
#'maxz=max(yest)#29.07#max(yest)
#'
#'dd=11
#'minx=min(Map[,1])
#'maxx=max(Map[,1])
#'miny=min(Map[,2])
#'maxy=max(Map[,2])
#'
#'se1=seq((minx+(maxx-minx)/(2*dd)),(maxx-(maxx-minx)/(dd*2)),(maxx-minx)/dd)
#'se2=seq((miny+(maxy-miny)/(2*dd)),(maxy-(maxy-miny)/(dd*2)),(maxy-miny)/dd)
#'
#'
#'graf=function(){
#'
#'for(rr in 1:nrow(yest)){
#'    vectMATRIZ=rep(maxz,dd*dd)
#'    vectMATRIZ[INDI]=yest[rr,]
#'    MATRIZ=matrix(vectMATRIZ,dd,dd)
#'    filled.contour(se1,se2,t(MATRIZ),nlevels=13,color = terrain.colors,zlim=c(infz,maxz),plot.axes = {
#'      lines(Map)},plot.title = title(main = paste(c('Year/Month/Season = ',YEAR[rr]),collapse=''),xlab = "Longitude", ylab = "Latitude"))
#'
#'  }
#'
#'}
#'
#'saveGIF(graf(),interval=0.5,movie.name="ModelA.gif")
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
