#' Defining the Model
#'
#' Give the complete formula and generate the expected kernel library.
#'
#' It processes data based on formula and label_names and creates a kernel
#' library according to the parameters given in Kern_par.
#'
#' @param formula A symbolic description of the model to be fitted.
#' @param label_names A character string indicating all the interior variables
#' included in each predictor.
#' @param data A dataframe to be fitted.
#' @param Kern_par A dataframe indicating the parameters of base kernels to be
#' created.
#' @return \item{Y}{Reponses of the dataframe.}
#'
#' \item{X1}{The first type of factor in the dataframe (could contains several
#' subfactors).}
#'
#' \item{X2}{The second type of factor in the dataframe (could contains several
#' subfactors).}
#'
#' \item{Kernlist}{The kernel library containing several kernels given by
#' user.}
#' @author Wenying Deng
#' @seealso method: \code{\link{kernelGenerate}}
#' @examples
#'
#'
#' ##Kern_par <- data.frame(method = c("rbf", "polynomial", "matern"),
#' ##Sigma = rep(0, 3), l = c(.5, 1, 1.5), p = 1:3)
#' ##Kern_par$method <- as.character(Kern_par$method)
#' ##defineModel(formula = Y ~ X1 + X2,
#' ##label_names = list(X1 = c("x1", "x2"), X2 = c("x3", "x4")),
#' ##data, Kern_par)
#'
#'
#' @export defineModel
defineModel <- function(formula, label_names, data, Kern_par){

  Y <- data[, as.character(attr(terms(formula), "variables"))[2]]
  re <- genericFormula(formula, label_names)
  generic_formula0 <- re$generic_formula
  len <- re$length_main

  X <- model.matrix(generic_formula0, data)[, -1]

  n <- nrow(X)
  Xm <- colMeans(X)
  p <- ncol(X)
  X <- X - rep(Xm, rep(n, p))
  Xscale <- drop(rep(1 / n, n) %*% X ^ 2) ^ .5
  X <- X / rep(Xscale, rep(n, p))

  X1 <- X[, c(1:length(label_names[[1]]))]
  X2 <- X[, c((length(label_names[[1]]) + 1):len)]

  Kernlist <- NULL
  for (d in 1:nrow(Kern_par))
    Kernlist <- c(Kernlist, kernelGenerate(Kern_par[d, ]$method,
                                           Kern_par[d, ]$Sigma,
                                           Kern_par[d, ]$l,
                                           Kern_par[d, ]$p))

  return(list(Y = Y, X1 = X1, X2 = X2, Kernlist = Kernlist))
}
