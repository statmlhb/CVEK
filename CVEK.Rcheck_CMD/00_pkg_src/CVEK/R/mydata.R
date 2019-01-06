#' Dataset as an example
#'
#' A dataset containing one column of fixed effects variable, two groups of
#' random effects variables, each with two columns of variables respectively. 
#' The variables are as follows:
#'
#' \itemize{
#'   \item Y. (vector of length n) Reponses of the dataframe.
#'   \item x1. (dataframe, n*1) Fixed effects variables in the dataframe.
#'   \item z1. (dataframe, n*1) The first column of the first group of random effects 
#'   variables in the dataframe.
#'   \item z2. (dataframe, n*1) The second column of the first group of random effects 
#'   variables in the dataframe.
#'   \item z3. (dataframe, n*1) The first column of the second group of random effects 
#'   variables in the dataframe.
#'   \item z4. (dataframe, n*1) The second column of the second group of random effects 
#'   variables in the dataframe.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name mydata
#' @usage data(mydata)
#' @format A data frame with 100 rows and 5 variables
NULL