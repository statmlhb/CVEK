#' Conducting Gaussian Process Regression
#' 
#' Conduct gaussian process regression based on the estimated ensemble kernel
#' matrix.
#' 
#' After obtaining the ensemble kernel matrix, we can calculate the outpur of
#' gaussian process regression.
#' 
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X (dataframe, n*p) Fixed effects variables in the dataframe (could
#' contains several subfactors).
#' @param Z1 (dataframe, n*q1) The first group of random effects variables in
#' the dataframe (could contains several subfactors).
#' @param Z2 (dataframe, n*q2) The second group of random effects variables in
#' the dataframe (could contains several subfactors).
#' @param kern_list (list of length K) A list of kernel functions given by
#' user.
#' @param mode (character) A character string indicating which tuning parameter
#' criteria is to be used.
#' @param strategy (character) A character string indicating which ensemble
#' strategy is to be used.
#' @param beta_exp (numeric/character) A numeric value specifying the parameter
#' when strategy = "exp" \code{\link{ensemble_exp}}.
#' @param lambda (numeric) A numeric string specifying the range of noise to be
#' chosen. The lower limit of lambda must be above 0.
#' @return \item{lambda}{(numeric) The selected tuning parameter based on the
#' estimated ensemble kernel matrix.}
#' 
#' \item{beta}{(matrix, p*1) Fixed effects estimator of the model.}
#' 
#' \item{alpha}{(vector of length n) Random effects estimator of the estimated
#' ensemble kernel matrix.}
#' 
#' \item{K}{(matrix, n*n) Estimated ensemble kernel matrix.}
#' 
#' \item{u_hat}{(vector of length K) A vector of weights of the kernels in the
#' library.}
#' 
#' \item{base_est}{(list of length 6) The detailed estimation results of K
#' kernels.}
#' @author Wenying Deng
#' @seealso strategy: \code{\link{ensemble}}
#' @examples
#' 
#' 
#' 
#' estimation(Y, X, Z1, Z2, kern_list, mode = "loocv", strategy = "stack",
#' beta_exp = 1, lambda = exp(seq(-10, 5)))
#' 
#' 
#' 
#' @export estimation
estimation <- function(Y, X, Z1, Z2, kern_list,
                       mode = "loocv", strategy = "stack", beta_exp = 1,
                       lambda = exp(seq(-10, 5))) {
  
  # base model estimate
  n <- length(Y)
  kern_size <- length(kern_list)
  base_est <- estimate_base(n, kern_size, Y, X, Z1, Z2, kern_list, 
                            mode, lambda)
  
  # ensemble estimate
  P_K_hat <- base_est$P_K_hat
  error_mat <- base_est$error_mat
  ens_res <- ensemble(n, kern_size, strategy, beta_exp, error_mat, P_K_hat)
  
  # assemble ensemble kernel matrix
  K_ens <- ensemble_kernel_matrix(ens_res$A_est)
  
  # final estimate
  lambda_ens <- tuning(Y, X, K_ens, mode, lambda)
  ens_est <- estimate_ridge(X = cbind(matrix(1, nrow = n, ncol = 1), X), 
                            K = K_ens, Y = Y, lambda = lambda_ens)
  list(lambda = lambda_ens, beta = ens_est$beta, 
       alpha = ens_est$alpha, K = K_ens, 
       u_hat = ens_res$u_hat, base_est = base_est)
}



#' Estimating Projection Matrices
#' 
#' Calculate the estiamted projection matrices for every kernels in the kernel
#' library.
#' 
#' For a given mode, this function return a list of projection matrices for
#' every kernels in the kernel library and a n*K matrix indicating errors.
#' 
#' @param n (integer) A numeric number specifying the number of observations.
#' @param kern_size (integer, equals to K) A numeric number specifying the
#' number of kernels in the kernel library.
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X (dataframe, n*p) Fixed effects variables in the dataframe (could
#' contains several subfactors).
#' @param Z1 (dataframe, n*q1) The first group of random effects variables in
#' the dataframe (could contains several subfactors).
#' @param Z2 (dataframe, n*q2) The second group of random effects variables in
#' the dataframe (could contains several subfactors).
#' @param kern_list (list of length K) A list of kernel functions given by
#' user.
#' @param mode (character) A character string indicating which tuning parameter
#' criteria is to be used.
#' @param lambda (numeric) A numeric string specifying the range of noise to be
#' chosen. The lower limit of lambda must be above 0.
#' @return \item{A_hat}{(list of length K) A list of projection matrices for
#' each kernel in the kernel library.}
#' 
#' \item{P_K_hat}{(list of length K) A list of projection matrices to kernel
#' space for each kernel in the kernel library.}
#' 
#' \item{beta_list}{(list of length K) A list of fixed effects estimator for
#' each kernel in the kernel library.}
#' 
#' \item{alpha_list}{(list of length K) A list of random effects estimator for
#' each kernel in the kernel library.}
#' 
#' \item{lambda}{(list of length K) A list of selected tuning parameters for
#' each kernel in the kernel library.}
#' 
#' \item{error_mat}{(matrix, n*K) A n\*K matrix indicating errors.}
#' @author Wenying Deng
#' @references Jeremiah Zhe Liu and Brent Coull. Robust Hypothesis Test for
#' Nonlinear Effect with Gaus- sian Processes. October 2017.
#' @examples
#' 
#' 
#' 
#' estimate_base(n = 100, kern_size = 3, Y, X, Z1, Z2, kern_list,
#' mode = "loocv", lambda = exp(seq(-10, 5)))
#' 
#' 
#' 
#' @export estimate_base
estimate_base <- function(n, kern_size, Y, X, Z1, Z2, kern_list, mode, lambda){
  A_hat <- list()
  P_K_hat <- list()
  beta_list <- list()
  alpha_list <- list()
  lambda_list <- list()
  error_mat <- matrix(0, nrow = n, ncol = kern_size)
  
  for (k in seq(kern_size)) {
    kern <- kern_list[[k]]
    K1_m <- kern(Z1, Z1)
    K2_m <- kern(Z2, Z2)
    K <- K1_m + K2_m
    K_scale <- tr(K)
    K <- K / K_scale
    
    if (length(lambda) != 0) {
      lambda0 <- tuning(Y, X, K, mode, lambda)
      estimate <- estimate_ridge(X = cbind(matrix(1, nrow = n, ncol = 1), X), 
                                 K = K, Y = Y, lambda = lambda0)
      A <- estimate$proj_matrix$total
      
      # produce loocv error matrix
      error_mat[, k] <- (diag(n) - A) %*% Y / (1 - diag(A))
      
      A_hat[[k]] <- A
      P_K_hat[[k]] <- estimate$proj_matrix$P_K0
      beta_list[[k]] <- estimate$beta
      alpha_list[[k]] <- estimate$alpha
      lambda_list[[k]] <- lambda0
    }
  }
  
  list(A_hat = A_hat, P_K_hat = P_K_hat,
       beta_list = beta_list, alpha_list = alpha_list,
       lambda_list = lambda_list, error_mat = error_mat)
}



#' Estimating a Single Model
#' 
#' Estimating projection matrices and parameter estimates for a single model.
#' 
#' For a single model, we can calculate the output of gaussian process
#' regression, the solution is given by \deqn{\hat{\beta}=[X^T(K+\lambda
#' I)^{-1}X]^{-1}X^T(K+\lambda I)^{-1}y} \deqn{\hat{\alpha}=(K+\lambda
#' I)^{-1}(y-\hat{\beta}X)}.
#' 
#' @param X (dataframe, n*p) Fixed effects variables in the dataframe (could
#' contains several subfactors).
#' @param K (matrix, n*n) Kernel matrix.
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param lambda (numeric) A numeric string specifying the range of noise to be
#' chosen. The lower limit of lambda must be above 0.
#' @return \item{lambda}{(numeric) The selected tuning parameter based on the
#' estimated ensemble kernel matrix.}
#' 
#' \item{beta}{(matrix, p*1) Fixed effects estimator of the model.}
#' 
#' \item{alpha}{(vector of length n) Random effects estimator of the estimated
#' ensemble kernel matrix.}
#' 
#' \item{proj_matrix}{(list of length 4) Estimated projection matrices of the
#' model.}
#' @author Wenying Deng
#' @examples
#' 
#' 
#' 
#' estimate_ridge(X = cbind(matrix(1, nrow = 100, ncol = 1), X), K_ens, Y, 
#' lambda = exp(seq(-10, 5)))
#' 
#' 
#' 
#' @export estimate_ridge
estimate_ridge <- function(X, K, Y, lambda){
  # standardize kernel matrix
  n <- length(Y)
  V_inv <- ginv(K + lambda * diag(n))
  B_mat <- ginv(t(X) %*% V_inv %*% X) %*% t(X) %*% V_inv
  
  # project matrices
  P_X <- X %*% B_mat  # projection to fixed-effect space
  P_K0 <- K %*% V_inv # projection to kernel space
  P_K <- P_K0 %*% (diag(n) - P_X) # residual projection to kernel space
  
  # parameter estimates
  beta <- B_mat %*% Y
  alpha <- V_inv %*% (diag(n) - P_X) %*% Y
  
  proj_matrix_list <- list(total = P_X + P_K, 
                           P_X = P_X, P_K = P_K, P_K0 = P_K0)
  list(beta = beta, alpha = alpha, 
       proj_matrix = proj_matrix_list)
}




#' Calculating Ensemble Kernel Matrix
#' 
#' Calculating ensemble kernel matrix and truncating those columns whose
#' eigenvalues are smaller than the given threshold.
#' 
#' After we obtain the ensemble projection matrix, we can calculate the
#' ensemble kernel matrix.
#' 
#' @param A_est (matrix) Ensemble projection matrix.
#' @param eig_thres (numeric) Threshold to truncate the kernel matrix.
#' @return \item{K_hat}{(matrix, n*n) Estimated ensemble kernel matrix.}
#' @author Wenying Deng
#' @examples
#' 
#' 
#' 
#' ensemble_kernel_matrix(A_hat[[2]], eig_thres = 1e-11)
#' 
#' 
#' 
#' @export ensemble_kernel_matrix
ensemble_kernel_matrix <- function(A_est, eig_thres = 1e-11){
  As <- svd(A_est)
  U <- As$u
  d <- As$d
  
  # produce spectral components for ensemble kernel matrix
  ensemble_dim <- sum(d > eig_thres)
  U_ens <- U[, 1:ensemble_dim]
  d_ens <- d[1:ensemble_dim]/(1 - d[1:ensemble_dim])
  
  # assemble ensemble matrix and return
  K_hat <- U_ens %*% diag(d_ens) %*% t(U_ens)
  K_hat / tr(K_hat)
}

