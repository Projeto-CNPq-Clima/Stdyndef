#' Summary Function for SpatialDeformationMCMC Output
#'
#' This function takes the output of the SpatialDeformationMCMC function and returns
#' a summary table with mean, 2.5% quantile, and 97.5% quantile for sigmak, phik, and pik.
#'
#' @param MCMC_output A list containing the matrices Msigmak, MPhi, and Mkappa, which are
#' the outputs from the SpatialDeformationMCMC function.
#' @return A data frame with summary statistics (mean, 2.5%, and 97.5%) for each parameter.
#' @export
summary_stdyndef <- function(MCMC_output) {


  # Extract parameters from MCMC output
  Msigmak <- MCMC_output$Msigmak
  MPhi <- MCMC_output$MPhi
  Mkappa <- MCMC_output$Mkappa

  # Calculate summary statistics
  summary_stats <- function(matrix) {
    mean_vals <- apply(matrix, 2, mean)
    lower_vals <- apply(matrix, 2, function(x) quantile(x, probs = 0.025))
    upper_vals <- apply(matrix, 2, function(x) quantile(x, probs = 0.975))

    data.frame(
      mean = mean_vals,
      `2.5%` = lower_vals,
      `97.5%` = upper_vals
    )
  }

  # Summarize each parameter
  sigmak_summary <- summary_stats(Msigmak)
  phik_summary <- summary_stats(MPhi)
  pik_summary <- summary_stats(Mkappa)

  # Combine summaries into a single data frame
  k_values <- 1:ncol(Msigmak)
  summary_table <- data.frame(
    k = k_values,
    `sigmak_mean` = sigmak_summary[,1],
    `sigmak_2.5%` = sigmak_summary[,2],
    `sigmak_97.5%` = sigmak_summary[,3],
    `phik_mean` = phik_summary[,1],
    `phik_2.5%` = phik_summary[,2],
    `phik_97.5%` = phik_summary[,3],
    `pik_mean` = pik_summary[,1],
    `pik_2.5%` = pik_summary[,2],
    `pik_97.5%` = pik_summary[,3]
  )

  return(summary_table)
}

# Exemplo de uso:
# result <- summary_SpatialDeformation(MCMC_output)
# print(result)
