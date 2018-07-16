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
#' @param Y Reponses of the dataframe.
#' @param X1 The first type of factor in the dataframe (could contains several
#' subfactors).
#' @param X2 The second type of factor in the dataframe (could contains several
#' subfactors).
#' @param Kernlist The kernel library containing several kernels given by user.
#' @param mode A character string indicating which tuning parameter criteria is
#' to be used.
#' @param strategy A character string indicating which ensemble strategy is to
#' be used.
#' @param beta A numeric value specifying the parameter when strategy = "exp".
#' @param lambda A numeric string specifying the range of noise to be chosen.
#' The lower limit of lambda must be above 0.
#' @return \item{lam}{The selected tuning parameter based on the estimated
#' ensemble kernel matrix.}
#'
#' \item{intercept}{Estimated bias of the model.}
#'
#' \item{alpha}{Estimated coefficients of the estimated ensemble kernel
#' matrix.}
#'
#' \item{K}{Estimated ensemble kernel matrix.}
#'
#' \item{u_hat}{A vector of weights of the kernels in the library.}
#' @author Wenying Deng
#' @seealso strategy: \code{\link{ensemble}}
#' @examples
#'
#'
#' ##estimation(Y, X1, X2, Kernlist, mode = "loocv", strategy = "erm",
#' ##beta = 1, lambda = exp(seq(-5, 5)))
#'
#'
#' @export estimation
estimation <- function(Y, X1, X2, Kernlist,
                       mode = "loocv", strategy = "erm", beta = 1,
                       lambda = exp(seq(-5, 5))){

  n <- length(Y)
  D <- length(Kernlist)
  out <- baseEstimate(n, D, Y, X1, X2, Kernlist, mode, lambda)
  A_hat <- out$A_hat
  error_mat <- out$error_mat

  out2 <- ensemble(n, D, strategy, beta, error_mat, A_hat)
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

  return(list(lam = lambda0, intercept = beta0,
              alpha = alpha, K = K_hat, u_hat = u_hat))
}
