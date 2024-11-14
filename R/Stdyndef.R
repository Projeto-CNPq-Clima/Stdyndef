#' Spatiotemporal Modeling with Dynamic Deformation for Nonstationary Covariance Structures
#'
#' The `Stdyndef` package provides tools for implementing and analyzing spatiotemporal models
#' with dynamic deformation, designed to handle nonstationary covariance structures. This innovative
#' modeling approach allows for dynamic variation in spatial correlation over time, enabling more
#' accurate and flexible analysis of complex spatiotemporal data.
#'
#' @docType package
#' @name Stdyndef
#'
#' @title Spatiotemporal Modeling with Dynamic Deformation for Nonstationary Covariance Structures
#'
#' @author
#' Fidel Ernesto Castro Morales \email{fidel.castro@ufrn.br} \cr
#' Marina S. Paez \email{marina@im.ufrj.br}
#'
#' @description
#' The `Stdyndef` package introduces a novel spatiotemporal modeling framework that incorporates
#' dynamic deformation to account for changes in the spatial correlation structure over time. By
#' leveraging state-space models, the package allows for smooth temporal deformation relative to the
#' original spatial region, making it suitable for data with evolving spatial relationships.
#'
#' @details
#' This package was developed to address the need for modeling monthly average temperature data in the
#' southern region of Brazil. The distinct geographic characteristics of this region, including plateaus,
#' mountain ranges, and proximity to the Atlantic Ocean, along with varying meteorological phenomena,
#' contribute to the need for a flexible and dynamic approach to spatial correlation modeling.
#'
#' The model parameters are estimated using a Bayesian framework, employing Markov Chain Monte Carlo
#' (MCMC) methods to approximate the posterior distribution of the parameters. This approach allows
#' users to:
#' - Model dynamic deformation of the spatial correlation structure.
#' - Analyze nonstationary spatiotemporal data.
#' - Perform Bayesian parameter estimation and uncertainty quantification.
#'
#' The package has been applied to 15 years of monthly average temperature data from the southern region
#' of Brazil. Results demonstrate significant improvements in temperature modeling using dynamic
#' deformation compared to static models, with the new approach capturing annual changes in correlation
#' structure.
#'
#' @examples
#' # Load the package and data
#' library(Stdyndef)
#' data(temperature)
#' data(FT)
#' data(MatFFT)
#' data(GT)
#' data(sites)
#' data(GAMA)
#'
#' # Fit the spatiotemporal model using MCMC
#' Mod <- SpatialDeformationMCMC(
#'   response = temperature,
#'   FT = FT,
#'   MatFFT = MatFFT,
#'   GT = GT,
#'   sites = sites,
#'   GAMA = GAMA,
#'   iteration = 1000,
#'   burnin = 500,
#'   jump = 10
#' )
#'
#' # Compute dynamic deformations
#' Mmd <- compute_Mmd(Mod, temperature, sites, Map, dd = 10)
#'
#' # Calculate the Deviance Information Criterion (DIC)
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
#' print(dic_results)
#'
#' @abstract
#' In this paper, we present an innovative spatiotemporal model that notably allows dynamic variation
#' in the spatial correlation structure over time through dynamic deformation. We propose that temporal
#' deformation occurs smoothly relative to the original region. To incorporate this idea, we employed
#' state-space models to model dynamic deformation. Generalizing the class of models that use spatial
#' deformation was driven by the need to model monthly average temperature data in the southern region
#' of Brazil. The distinctive traits of this region, characterized by plateaus and mountain ranges and
#' close proximity to the Atlantic Ocean, give rise to geographic diversity. This diversity, in addition
#' to different meteorological phenomena over time, may influence the spatial correlation function.
#'
#' The model parameters were estimated using a Bayesian approach, requiring the use of MCMC methods to
#' approximate the posterior distribution of the parameters. The model was applied to 15 years of monthly
#' average temperature data from the southern region of Brazil. The primary result of this analysis
#' revealed a significant improvement in temperature modeling when using the proposed model compared to
#' that using versions that employ static deformation. We believe this improvement is due to the ability
#' of the new model to allow annual changes in the correlation structure.
