#' Conducting Gaussian Process Regression
#'
#' Conduct gaussian process regression based on the estimated ensemble kernel
#' matrix.
#'
#' After obtaining the ensemble kernel matrix, we can calculate the outpur of
#' gaussian process regression, the solution is given by
#' \deqn{\hat{\beta}=[1^T(K+\lambda I)^{-1}1]^{-1}1^T(K+\lambda I)^{-1}y}
#' \deqn{\hat{\alpha}=(K+\lambda I)^{-1}(y-\hat{\beta}1)} where
#' \eqn{\beta=intercept}.
#'
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X1 (dataframe, n*p1) The first type of factor in the dataframe (could 
#' contains several subfactors).
#' @param X2 (dataframe, n*p2) The second type of factor in the dataframe (could 
#' contains several subfactors).
#' @param kern_list (list of length K) A list of kernel functions given by user.
#' @param mode (character) A character string indicating which tuning 
#' parameter criteria is to be used.
#' @param strategy (character) A character string indicating which ensemble 
#' strategy is to be used.
#' @param beta (numeric/character) A numeric value specifying the parameter 
#' when strategy = "exp" \code{\link{ensemble_exp}}.
#' @param lambda (numeric) A numeric string specifying the range of 
#' noise to be chosen. The lower limit of lambda must be above 0.
#' @return \item{lam}{(numeric) The selected tuning parameter based on the 
#' estimated ensemble kernel matrix.}
#'
#' \item{intercept}{(numeric) Estimated bias of the model.}
#'
#' \item{alpha}{(vector of length n) Estimated coefficients of the estimated 
#' ensemble kernel matrix.}
#'
#' \item{K}{(matrix, n*n) Estimated ensemble kernel matrix.}
#'
#' \item{u_hat}{(vector of length K) A vector of weights of the kernels in 
#' the library.}
#' @author Wenying Deng
#' @seealso strategy: \code{\link{ensemble}}
#' @examples
#'
#'
#' estimation(Y, X1, X2, kern_list, mode = "loocv", strategy = "erm",
#' beta = 1, lambda = exp(seq(-5, 5)))
#'
#'
#' @export estimation
estimation <- function(Y, X1, X2, kern_list,
                       mode = "loocv", strategy = "erm", beta = 1,
                       lambda = exp(seq(-5, 5))) {
  
  n <- length(Y)
  kern_size <- length(kern_list)
  out <- estimate_base(n, kern_size, Y, X1, X2, kern_list, mode, lambda)
  A_hat <- out$A_hat
  error_mat <- out$error_mat
  
  out2 <- ensemble(n, kern_size, strategy, beta, error_mat, A_hat)
  A_est <- out2$A_est
  u_hat <- out2$u_hat
  
  As <- svd(A_est)
  K_hat <- As$u %*% diag(As$d / (1 - As$d)) %*% t(As$u)
  
  lambda0 <- tuning(Y, K_hat, mode, lambda)
  K1 <- cbind(1, K_hat)
  K2 <- cbind(0, rbind(0, K_hat))
  
  theta <- ginv(lambda0 * K2 + t(K1) %*% K1) %*% t(K1) %*% Y
  beta0 <- theta[1]
  alpha <- theta[-1]
  
  list(lam = lambda0, intercept = beta0, 
       alpha = alpha, K = K_hat, u_hat = u_hat)
}

#' Estimating Projection Matrices
#' 
#' Calculate the estiamted projection matrices for every kernels in the kernel
#' library.
#' 
#' For a given mode, this function return a list of projection matrices for
#' every kernels in the kernel library and a n*kern_size matrix indicating
#' errors.
#' 
#' @param n (integer) A numeric number specifying the number of observations.
#' @param kern_size (integer, =K) A numeric number specifying the number of 
#' kernels in the kernel library.
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X1 (dataframe, n*p1) The first type of factor in the dataframe (could 
#' contains several subfactors).
#' @param X2 (dataframe, n*p2) The second type of factor in the dataframe (could 
#' contains several subfactors).
#' @param kern_list (list of length K) A list of kernel functions given by user.
#' @param mode (character) A character string indicating which tuning parameter 
#' criteria is to be used.
#' @param lambda (numeric) A numeric string specifying the range of noise 
#' to be chosen. The lower limit of lambda must be above 0.
#' @return \item{A_hat}{(list of length K) A list of projection matrices for 
#' every kernels in the kernel library.}
#' 
#' \item{error_mat}{(matrix, n*K) A n\*kern_size matrix indicating errors.}
#' @author Wenying Deng
#' @references Jeremiah Zhe Liu and Brent Coull. Robust Hypothesis Test for
#' Nonlinear Effect with Gaus- sian Processes. October 2017.
#' @examples
#' 
#' 
#' estimate_base(n = 100, kern_size = 3, Y, X1, X2, kern_list,
#' mode = "loocv", lambda = exp(seq(-5, 5)))
#' 
#' 
#' @export estimate_base
estimate_base <- function(n, kern_size, Y, X1, X2, 
                          kern_list, mode, lambda) {
  
  A_hat <- list()
  error_mat <- matrix(0, nrow = n, ncol = kern_size)
  
  for (d in seq(kern_size)) {
    kern <- kern_list[[d]]
    K1_m <- kern(X1, X1)
    K2_m <- kern(X2, X2)
    if (tr(K1_m) > 0 & tr(K2_m) > 0) {
      K1_m <- K1_m / tr(K1_m)
      K2_m <- K2_m / tr(K2_m)
    }
    K <- K1_m + K2_m
    if (length(lambda) != 1) {
      lambda0 <- tuning(Y, K, mode, lambda)
      K1 <- cbind(1, K)
      K2 <- cbind(0, rbind(0, K))
      theta <- ginv(lambda0 * K2 + t(K1) %*% K1) %*% t(K1) %*% Y
      beta0 <- theta[1]
      M <- K %*% ginv(K + lambda0 * diag(n))
      error_mat[, d] <- (diag(n) - M) %*% (Y - beta0) / diag(diag(n) - M)
      A_hat[[d]] <- M
    }
  }
  
  list(A_hat = A_hat, error_mat = error_mat)
}
