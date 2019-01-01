#' Defining the Model
#' 
#' Give the complete formula and generate the expected kernel library.
#' 
#' It processes data based on formula and label_names and creates a kernel
#' library according to the parameters given in kern_par.
#' 
#' * label_names: for two groups of random effects with sizes q1 and q2
#'   respectively, label_names contains two elements. The length of the first
#'   element is q1, indicating the names of q1 interiors variables, and the
#'   length of second one is q2, indicating the names of q2 interiors variables.
#' 
#' * data: for a data with n observations and p+q variables (with sub-groups of
#'   sizes (q1, q2), q=q1+q2), the dimension of dataframe is n*(p+q). All entries
#'   should be numeric and the column name of response is "Y", while the column
#'   names of q variables are the ones from label_names.
#' 
#' * kern_par: for a library of K kernels, the dimension of this dataframe is
#'   K*3. Each row represents a kernel. The first column is method, with entries
#'   of character class. The second and the third are l and p respectively, both
#'   with entries of numeric class.
#' 
#' @param formula (formula) A symbolic description of the model to be fitted.
#' @param data (dataframe, n*(p+q1+q2)) A dataframe to be fitted. See Details.
#' @param kern_par (dataframe, K*3) A dataframe indicating the parameters of
#' base kernels to fit random effects. See Details.
#' @param fixed_num (integer) A numeric number specifying the dimension of
#' fixed effects.
#' @param label_names (list) A character string indicating all the interior
#' variables included in each group of random effect. See Details.
#' @return \item{Y}{(vector of length n) Reponses of the dataframe.}
#' 
#' \item{X}{(dataframe, n*p) Fixed effects variables in the dataframe (could
#' contains several subfactors).}
#' 
#' \item{Z1}{(dataframe, n*q1) The first group of random effects variables in
#' the dataframe (could contains several subfactors).}
#' 
#' \item{Z2}{(dataframe, n*q2) The second group of random effects variables in
#' the dataframe (could contains several subfactors).}
#' 
#' \item{kern_list}{(list of length K) A list of kernel functions given by
#' user.}
#' @author Wenying Deng
#' @seealso method: \code{\link{generate_kernel}}
#' @examples
#' 
#' 
#' 
#' kern_par <- data.frame(method = c("rbf", "polynomial", "matern"), 
#' l = c(.5, 1, 1.5), d = 1:3)
#' kern_par$method <- as.character(kern_par$method)
#' define_model(formula = Y ~ X + Z1 + Z2, data = mydata, kern_par, 
#' fixed_num = 1, label_names = list(Z1 = c("z1", "z2"), Z2 = c("z3", "z4")))
#' 
#' 
#' 
#' @export define_model
define_model <- function(formula, data, kern_par = NULL, 
                         fixed_num = 1, label_names = NULL) {
  
  if ((fixed_num == 0) & is.null(label_names)) {
    stop("fixed effect and random effect can not be null simultaneously!")
  }
  Y <- data[, as.character(attr(terms(formula), "variables"))[2]]
  n <- length(Y)
  kern_list <- list()
  X <- NULL
  Z1 <- NULL
  Z2 <- NULL
  kern_list <- NULL
  if (!is.null(label_names)) {
    re <- generate_formula(formula, label_names)
    generic_formula0 <- re$generic_formula
    len <- re$length_main
    Z <- model.matrix(generic_formula0, data)[, -1]
    Z <- standardize(Z)
    Z1 <- Z[, c(1:length(label_names[[1]]))]
    Z2 <- Z[, c((length(label_names[[1]]) + 1):len)]
    
    for (k in 1:nrow(kern_par)) {
      kern_list[[k]] <- generate_kernel(kern_par[k, ]$method,
                                        kern_par[k, ]$l,
                                        kern_par[k, ]$d)
    }
  }
  if (fixed_num > 0) {
    X <- as.matrix(data[, 2:(1 + fixed_num)], nrow = n)
    X <- standardize(X)
  }
  
  list(Y = Y, X = X, Z1 = Z1, Z2 = Z2, kern_list = kern_list)
}
