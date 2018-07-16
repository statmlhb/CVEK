#' Estimating Ensemble Kernel Matrices
#' 
#' Give a list of estimated kernel matrices and their weights.
#' 
#' There are three ensemble strategies available here:
#' 
#' \bold{Empirical Risk Minimization}
#' 
#' After obtaining the estimated errors \eqn{\{\hat{\epsilon}_d\}_{d=1}^D}, we
#' estimate the ensemble weights \eqn{u=\{u_d\}_{d=1}^D} such that it minimizes
#' the overall error \deqn{\hat{u}=\underset{u \in \Delta}{argmin}\parallel
#' \sum_{d=1}^Du_d\hat{\epsilon}_d\parallel^2 \quad where\; \Delta=\{u | u \geq
#' 0, \parallel u \parallel_1=1\}} Then produce the final ensemble prediction:
#' \deqn{\hat{h}=\sum_{d=1}^D \hat{u}_d h_d=\sum_{d=1}^D \hat{u}_d
#' A_{d,\hat{\lambda}_d}y=\hat{A}y} where \eqn{\hat{A}=\sum_{d=1}^D \hat{u}_d
#' A_{d,\hat{\lambda}_d}} is the ensemble matrix.
#' 
#' \bold{Simple Averaging}
#' 
#' Motivated by existing literature in omnibus kernel, we propose another way
#' to obtain the ensemble matrix by simply choosing unsupervised weights
#' \eqn{u_d=1/D} for \eqn{d=1,2,...D}.
#' 
#' \bold{Exponential Weighting}
#' 
#' Additionally, another scholar gives a new strategy to calculate weights
#' based on the estimated errors \eqn{\{\hat{\epsilon}_d\}_{d=1}^D}.
#' \deqn{u_d(\beta)=\frac{exp(-\parallel \hat{\epsilon}_d
#' \parallel_2^2/\beta)}{\sum_{d=1}^Dexp(-\parallel \hat{\epsilon}_d
#' \parallel_2^2/\beta)}}
#' 
#' @param n A numeric number specifying the number of observations.
#' @param D A numeric number specifying the number of kernels in the kernel
#' library.
#' @param strategy A character string indicating which ensemble strategy is to
#' be used.
#' @param beta A numeric value specifying the parameter when strategy = "exp".
#' @param error_mat A n*D matrix indicating errors.
#' @param A_hat A list of projection matrices for every kernels in the kernel
#' library.
#' @return \item{A_est}{A list of estimated kernel matrices.}
#' 
#' \item{u_hat}{A vector of weights of the kernels in the library.}
#' @author Wenying Deng
#' @seealso mode: \code{\link{tuning}}
#' @references Jeremiah Zhe Liu and Brent Coull. Robust Hypothesis Test for
#' Nonlinear Effect with Gaus- sian Processes. October 2017.
#' 
#' Xiang Zhan, Anna Plantinga, Ni Zhao, and Michael C. Wu. A fast small-sample
#' kernel inde- pendence test for microbiome community-level association
#' analysis. December 2017.
#' 
#' Arnak S. Dalalyan and Alexandre B. Tsybakov. Aggregation by Exponential
#' Weighting and Sharp Oracle Inequalities. In Learning Theory, Lecture Notes
#' in Computer Science, pages 97– 111. Springer, Berlin, Heidelberg, June 2007.
#' @examples
#' 
#' 
#' ##ensemble(n = 50, D = 6, strategy = "erm", beta = 1, error_mat, A_hat)
#' 
#' 
#' @export ensemble
ensemble <-
  function(n, D, strategy, beta, error_mat, A_hat){

    if (strategy == "erm"){

      A <- error_mat
      B <- rep(0, n)
      E <- rep(1, D)
      F <- 1
      G <- diag(D)
      H <- rep(0, D)
      u_hat <- lsei(A, B, E = E, F = F, G = G, H = H)$X
      A_est <- u_hat[1] * A_hat[[1]]
      if(D != 1)
        for(d in 2:D)
          A_est <- A_est + u_hat[d] * A_hat[[d]]
    }
    else if (strategy == "average"){

	  u_hat <- rep(1 / D, D)
      A_est <- (1 / D) * A_hat[[1]]
      if(D != 1)
        for(d in 2:D)
          A_est <- A_est + (1 / D) * A_hat[[d]]
    }
    else if (strategy == "exp"){

      A <- error_mat
      u_hat <- apply(A, 2, function(x) exp(sum(-x ^ 2 / beta)))
      u_hat <- u_hat / sum(u_hat)
      A_est <- u_hat[1] * A_hat[[1]]
      if(D != 1)
        for(d in 2:D)
          A_est <- A_est + u_hat[d] * A_hat[[d]]
    }
    else
      stop("strategy must be erm, average or exp!")

    return(list(A_est = A_est, u_hat = u_hat))
  }








#' Estimating Projection Matrices
#' 
#' Calculate the estiamted projection matrices for every kernels in the kernel
#' library.
#' 
#' For a given mode, this function return a list of projection matrices for
#' every kernels in the kernel library and a size*magn matrix indicating
#' errors.
#' 
#' @param size A numeric number specifying the number of observations.
#' @param magn A numeric number specifying the number of kernels in the kernel
#' library.
#' @param Y Reponses of the dataframe.
#' @param X1 The first type of factor in the dataframe (could contains several
#' subfactors).
#' @param X2 The second type of factor in the dataframe (could contains several
#' subfactors).
#' @param Kernlist The kernel library containing several kernels given by user.
#' @param mode A character string indicating which tuning parameter criteria is
#' to be used.
#' @param lambda A numeric string specifying the range of noise to be chosen.
#' The lower limit of lambda must be above 0.
#' @return \item{A_hat}{A list of projection matrices for every kernels in the
#' kernel library.}
#' 
#' \item{error_mat}{A size*magn matrix indicating errors.}
#' @author Wenying Deng
#' @references Jeremiah Zhe Liu and Brent Coull. Robust Hypothesis Test for
#' Nonlinear Effect with Gaus- sian Processes. October 2017.
#' @examples
#' 
#' 
#' ##baseEstimate(size = 50, magn = 3, Y, X1, X2, Kernlist = NULL,
#' ##mode = "loocv", lambda = exp(seq(-5, 5)))
#' 
#' 
#' @export baseEstimate
baseEstimate <- function(size, magn, Y, X1, X2, Kernlist, mode, lambda){

  A_hat <- list()
  error_mat <- matrix(0, nrow = size, ncol = magn)

  for (d in seq(magn)){
    Kern <- Kernlist[[d]]
    K1_m <- Kern(X1, X1)
    K2_m <- Kern(X2, X2)
    if(tr(K1_m) > 0 & tr(K2_m) > 0){
      K1_m <- K1_m / tr(K1_m)
      K2_m <- K2_m / tr(K2_m)
    }
    K <- K1_m + K2_m
    if (length(lambda) != 1){
      lambda0 <- tuning(Y, K, mode, lambda)
      K1 <- cbind(1, K)
      K2 <- cbind(0, rbind(0, K))
      theta <- ginv(lambda0 * K2 + t(K1) %*% K1) %*% t(K1) %*% Y
      beta0 <- theta[1]
      M <- K %*% ginv(K + lambda0 * diag(size))
      error_mat[, d] <- (diag(size) - M) %*% (Y - beta0) / diag(diag(size) - M)
      A_hat[[d]] <- M
    }
  }

  return(list(A_hat = A_hat, error_mat = error_mat))
}
