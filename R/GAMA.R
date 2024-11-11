#' Partition Matrix for Time Indices in State-Space Model
#'
#' The `GAMA` matrix defines a partition of time indices for use in the state-space model. With dimensions
#' `q x r` such that `q * r = T`, this matrix organizes the time index set `AT = {1, ..., T}` into subsets.
#' Each row of `GAMA` contains the indices that form a subset of `AT`, ensuring that each subset is unique
#' and non-overlapping, and that the union of all subsets covers the entire set `AT`.
#'
#' @format A matrix with `q` rows and `r` columns:
#' \describe{
#'   \item{q}{Number of subsets of `AT`.}
#'   \item{r}{Number of time indices within each subset, such that `q * r = T`.}
#' }
#'
#' @details
#' The `GAMA` matrix is structured to create a partition of the time index set `AT = {1, ..., T}`. Each row
#' `i` of `GAMA` corresponds to a subset `GAMA_i = AT[GAMA[i,]]` for `i = 1, ..., q`. This partitioning is
#' such that:
#' - Each subset `GAMA_i` is disjoint from all other subsets, i.e., `GAMA_i ∩ GAMA_j` is empty for all `i ≠ j`.
#' - The union of all subsets `GAMA_i` equals `AT`, ensuring complete coverage of the time indices.
#'
#' This partitioning is useful for spatiotemporal models where temporal indices are grouped for analysis
#' across different time segments or batches.
#'
#' @examples
#' # Load the GAMA matrix
#' data(GAMA)
#'
#' # Display the structure of GAMA
#' str(GAMA)
#'
"GAMA"
