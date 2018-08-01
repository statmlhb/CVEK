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
#' @param n (integer) A numeric number specifying the number of observations.
#' @param kern_size (integer, =K) A numeric number specifying the number of 
#' kernels in the kernel library.
#' @param strategy (character) A character string indicating which ensemble 
#' strategy is to be used.
#' @param beta (numeric/character) A numeric value specifying the parameter when 
#' strategy = "exp" \code{\link{ensemble_exp}}.
#' @param error_mat (matrix, n*K) A n\*kern_size matrix indicating errors.
#' @param A_hat (list of length K) A list of projection matrices for 
#' every kernels in the kernel library.
#' @return \item{A_est}{(matrix, n*n) A list of estimated kernel matrices.}
#' 
#' \item{u_hat}{(vector of length K) A vector of weights of the kernels 
#' in the library.}
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
#' ensemble(n = 100, kern_size = 3, strategy = "erm", beta = 1, 
#' error_mat, A_hat)
#' 
#' 
#' @export ensemble

ensemble <-
  function(n, kern_size, strategy, beta, error_mat, A_hat) {
    
    strategy <- match.arg(strategy, c("avg", "exp", "erm"))
    func_name <- paste0("ensemble_", strategy)
    do.call(func_name, list(n = n, kern_size = kern_size, beta = beta,
                            error_mat = error_mat, A_hat = A_hat))
  }



#' Estimating Ensemble Kernel Matrices Using ERM
#' 
#' Give a list of estimated kernel matrices and their weights using 
#' empirical risk minimization.
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
#' @param n (integer) A numeric number specifying the number of observations.
#' @param kern_size (integer, =K) A numeric number specifying the number of 
#' kernels in the kernel library.
#' @param beta (numeric/character) A numeric value specifying the parameter when 
#' strategy = "exp" \code{\link{ensemble_exp}}.
#' @param error_mat (matrix, n*K) A n\*kern_size matrix indicating errors.
#' @param A_hat (list of length K) A list of projection matrices for 
#' every kernels in the kernel library.
#' @return \item{A_est}{(matrix, n*n) A list of estimated kernel matrices.}
#' 
#' \item{u_hat}{(vector of length K) A vector of weights of the kernels 
#' in the library.}
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
#' 
#' @export ensemble_erm
ensemble_erm <- 
  function(n, kern_size, beta, error_mat, A_hat) {
    A <- error_mat
    B <- rep(0, n)
    E <- rep(1, kern_size)
    F <- 1
    G <- diag(kern_size)
    H <- rep(0, kern_size)
    u_hat <- lsei(A, B, E = E, F = F, G = G, H = H)$X
    A_est <- u_hat[1] * A_hat[[1]]
    if (kern_size != 1) {
      for (d in 2:kern_size) {
        A_est <- A_est + u_hat[d] * A_hat[[d]]
      }
    }
    
    list(A_est = A_est, u_hat = u_hat)
  }


#' Estimating Ensemble Kernel Matrices Using AVG
#' 
#' Give a list of estimated kernel matrices and their weights using 
#' simple averaging.
#' 
#' \bold{Simple Averaging}
#' 
#' Motivated by existing literature in omnibus kernel, we propose another way
#' to obtain the ensemble matrix by simply choosing unsupervised weights
#' \eqn{u_d=1/D} for \eqn{d=1,2,...D}.
#' 
#' @param n (integer) A numeric number specifying the number of observations.
#' @param kern_size (integer, =K) A numeric number specifying the number of 
#' kernels in the kernel library.
#' @param beta (numeric/character) A numeric value specifying the parameter when 
#' strategy = "exp" \code{\link{ensemble_exp}}.
#' @param error_mat (matrix, n*K) A n\*kern_size matrix indicating errors.
#' @param A_hat (list of length K) A list of projection matrices for 
#' every kernels in the kernel library.
#' @return \item{A_est}{(matrix, n*n) A list of estimated kernel matrices.}
#' 
#' \item{u_hat}{(vector of length K) A vector of weights of the kernels 
#' in the library.}
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
#' 
#' @export ensemble_avg
ensemble_avg <- 
  function(n, kern_size, beta, error_mat, A_hat) {
    
    u_hat <- rep(1 / kern_size, kern_size)
    A_est <- (1 / kern_size) * A_hat[[1]]
    if (kern_size != 1) {
      for (d in 2:kern_size) {
        A_est <- A_est + (1 / kern_size) * A_hat[[d]]
      }
    }
    
    list(A_est = A_est, u_hat = u_hat)
  }



#' Estimating Ensemble Kernel Matrices Using EXP
#' 
#' Give a list of estimated kernel matrices and their weights using
#' exponential weighting.
#' 
#' \bold{Exponential Weighting}
#' 
#' Additionally, another scholar gives a new strategy to calculate weights
#' based on the estimated errors \eqn{\{\hat{\epsilon}_d\}_{d=1}^D}.
#' \deqn{u_d(\beta)=\frac{exp(-\parallel \hat{\epsilon}_d
#' \parallel_2^2/\beta)}{\sum_{d=1}^Dexp(-\parallel \hat{\epsilon}_d
#' \parallel_2^2/\beta)}}
#' 
#' \bold{beta}
#' 
#' The value of beta can be "min"=\eqn{min\{RSS\}_{d=1}^D/10},
#' "med"=\eqn{median\{RSS\}_{d=1}^D}, "max"=\eqn{max\{RSS\}_{d=1}^D*2}
#' and any other positive numeric number, where \eqn{\{RSS\} _{d=1}^D}
#' are the set of residual sum of squares of \eqn{D} base kernels.
#' 
#' @param n (integer) A numeric number specifying the number of observations.
#' @param kern_size (integer, =K) A numeric number specifying the number of 
#' kernels in the kernel library.
#' @param beta (numeric/character) A numeric value specifying the parameter 
#' when strategy = "exp". See Details.
#' @param error_mat (matrix, n*K) A n\*kern_size matrix indicating errors.
#' @param A_hat (list of length K) A list of projection matrices for 
#' every kernels in the kernel library.
#' @return \item{A_est}{(matrix, n*n) A list of estimated kernel matrices.}
#' 
#' \item{u_hat}{(vector of length K) A vector of weights of the kernels 
#' in the library.}
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
#' 
#' @export ensemble_exp
ensemble_exp <- 
  function(n, kern_size, beta, error_mat, A_hat) {
    
    A <- error_mat
    
    if (beta == "med") {
      beta <- median(apply(A, 2, function(x) sum(x ^ 2)))
    } else if (beta == "min") {
      beta <- min(apply(A, 2, function(x) sum(x ^ 2))) / 10
    } else if (beta == "max") {
      beta <- max(apply(A, 2, function(x) sum(x ^ 2))) * 2
    }
    
    u_hat <- apply(A, 2, function(x) {
      exp(sum(-x ^ 2 / beta))
    })
    u_hat <- u_hat / sum(u_hat)
    A_est <- u_hat[1] * A_hat[[1]]
    if (kern_size != 1) {
      for (d in 2:kern_size) {
        A_est <- A_est + u_hat[d] * A_hat[[d]]
      }
    }
    
    list(A_est = A_est, u_hat = u_hat)
  }
