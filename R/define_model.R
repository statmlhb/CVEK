#' Defining the Model
#'
#' Give the complete formula and generate the expected kernel library.
#'
#' It processes data based on formula and label_names and creates a kernel
#' library according to the parameters given in kern_par.
#' 
#' * label_names: for two subgroups with sizes p1 and p2 respectively,
#'   label_names contains two elements. The length of the first element is p1, 
#'   indicating the names of p1 interiors variables, and the length of second 
#'   one is p2, indicating the names of p2 interiors variables.
#'   
#' * data: for a data with n observations and P=p1+p2 variables (with sub-groups 
#'   of sizes (p1, p2)), the dimension of dataframe is n*P. All entries should
#'   be numeric and the column name of response is "Y", while the column 
#'   names of P variables are the ones from label_names.
#'   
#' * kern_par: for a library of K kernels, the dimension of this dataframe is
#'   K*4. Each row represents a kernel. The first column is method, with entries
#'   of character class. The second is Sigma, with entries of matrix class, 
#'   indicating the covariance matrix for neural network kernel (default=0).
#'   The third and the fourth are l and p respectively, both with entries of 
#'   numeric class.
#'
#' @param formula (formula) A symbolic description of the model to be fitted.
#' @param label_names (list) A character string indicating all the interior 
#' variables included in each predictor. See Details.
#' @param data (dataframe, n*P) A dataframe to be fitted. See Details.
#' @param kern_par (dataframe, K*4) A dataframe indicating the parameters 
#' of base kernels to be created. See Details.
#' @return \item{Y}{(vector of length n) Reponses of the dataframe.}
#'
#' \item{X1}{(dataframe, n*p1) The first type of factor in the dataframe (could 
#' contains several subfactors).}
#'
#' \item{X2}{(dataframe, n*p2) The second type of factor in the dataframe (could 
#' contains several subfactors).}
#'
#' \item{kern_list}{(list of length K) A list of kernel functions given by user.}
#' @author Wenying Deng
#' @seealso method: \code{\link{generate_kernel}}
#' @examples
#'
#'
#' kern_par <- data.frame(method = c("rbf", "polynomial", "matern"), 
#' Sigma = rep(0, 3), l = c(.5, 1, 1.5), p = 1:3)
#' kern_par$method <- as.character(kern_par$method)
#' define_model(formula = Y ~ X1 + X2,
#' label_names = list(X1 = c("x1", "x2"), X2 = c("x3", "x4")),
#' data = dora, kern_par)
#'
#'
#' @export define_model
define_model <- function(formula, label_names, data, kern_par) {
  
  Y <- data[, as.character(attr(terms(formula), "variables"))[2]]
  re <- generate_formula(formula, label_names)
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
  
  kern_list <- list()
  for (d in 1:nrow(kern_par)) {
    kern_list[[d]] <- generate_kernel(kern_par[d, ]$method,
                                      kern_par[d, ]$Sigma,
                                      kern_par[d, ]$l,
                                      kern_par[d, ]$p)
  }
  
  list(Y = Y, X1 = X1, X2 = X2, kern_list = kern_list)
}
