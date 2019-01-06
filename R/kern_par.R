#' Dataframe as an example for kernel parameters
#'
#' A dataframe containing the information from user to define base kernels.
#'
#' \itemize{
#'   \item method. (character) A character string indicating which kernel 
#' is to be computed.
#'   \item l. (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#'   \item d. (integer) For polynomial, d is the power; for matern, v = d + 1 / 2; for
#' rational, alpha = d.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name kern_par
#' @usage data(kern_par)
#' @format A data frame with 3 rows and 3 variables
NULL