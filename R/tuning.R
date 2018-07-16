#' Calculating Tuning Parameters
#' 
#' Calculate tuning parameters based on given criteria.
#' 
#' There are four tuning parameter selections here:
#' 
#' \bold{leave-one-out Cross Validation}
#' 
#' \deqn{\lambda_{n-CV}=\underset{\lambda \in
#' \Lambda}{argmin}\;\Big\{log\;y^{\star
#' T}[I-diag(A_\lambda)-\frac{1}{n}I]^{-1}(I-A_\lambda)^2[I-diag(A_\lambda)-\frac{1}{n}I]^{-1}y^\star
#' \Big\}}
#' 
#' \bold{Akaike Information Criteria}
#' 
#' \deqn{\lambda_{AICc}=\underset{\lambda \in \Lambda}{argmin}\Big\{log\;
#' y^{\star
#' T}(I-A_\lambda)^2y^\star+\frac{2[tr(A_\lambda)+2]}{n-tr(A_\lambda)-3}\Big\}}
#' 
#' \bold{Generalized Cross Validation}
#' 
#' \deqn{\lambda_{GCVc}=\underset{\lambda \in \Lambda}{argmin}\Big\{log\;
#' y^{\star
#' T}(I-A_\lambda)^2y^\star-2log[1-\frac{tr(A_\lambda)}{n}-\frac{2}{n}]_+\Big\}}
#' 
#' \bold{Generalized Maximum Profile Marginal Likelihood}
#' 
#' \deqn{\lambda_{GMPML}=\underset{\lambda \in \Lambda}{argmin}\Big\{log\;
#' y^{\star T}(I-A_\lambda)y^\star-\frac{1}{n-1}log \mid I-A_\lambda \mid
#' \Big\}}
#' 
#' @param Y Reponses of the dataframe.
#' @param K_mat Estimated ensemble kernel matrix.
#' @param mode A character string indicating which tuning parameter criteria is
#' to be used.
#' @param lambda A numeric string specifying the range of noise to be chosen.
#' The lower limit of lambda must be above 0.
#' @return \item{lambda0}{The estimated tuning parameter.}
#' @author Wenying Deng
#' @references Philip S. Boonstra, Bhramar Mukherjee, and Jeremy M. G. Taylor.
#' A Small-Sample Choice of the Tuning Parameter in Ridge Regression. July
#' 2015.
#' 
#' Trevor Hastie, Robert Tibshirani, and Jerome Friedman. The Elements of
#' Statistical Learning: Data Mining, Inference, and Prediction, Second
#' Edition. Springer Series in Statistics. Springer- Verlag, New York, 2
#' edition, 2009.
#' 
#' Hirotogu Akaike. Information Theory and an Extension of the Maximum
#' Likelihood Princi- ple. In Selected Papers of Hirotugu Akaike, Springer
#' Series in Statistics, pages 199–213. Springer, New York, NY, 1998.
#' 
#' Clifford M. Hurvich and Chih-Ling Tsai. Regression and time series model
#' selection in small samples. June 1989.
#' 
#' Hurvich Clifford M., Simonoff Jeffrey S., and Tsai Chih-Ling. Smoothing
#' parameter selection in nonparametric regression using an improved Akaike
#' information criterion. January 2002.
#' @examples
#' 
#' 
#' ##tuning(Y, K_mat = K, mode = "loocv", lambda = exp(seq(-5, 5)))
#' 
#' 
#' @export tuning
tuning <-
  function(Y, K_mat, mode, lambda){

    n <- nrow(K_mat)

    # estimation
    if (mode == "loocv"){
      CV <- sapply(lambda, function(k){
        A <- K_mat %*% ginv(K_mat + k * diag(n))
        sum(((diag(n) - A) %*% Y / diag(diag(n) - A)) ^ 2)
      })
    }
    else if (mode == "AICc"){
      CV <- sapply(lambda, function(k){
        A <- K_mat %*% ginv(K_mat + k * diag(n))
        log(t(Y) %*% (diag(n) - A) %*% (diag(n) - A) %*% Y) +
          2 * (tr(A) + 2) / (n - tr(A) - 3)
      })
    }
    else if (mode == "GCVc"){
      CV <- sapply(lambda, function(k){
        A <- K_mat %*% ginv(K_mat + k * diag(n))
        log(t(Y) %*% (diag(n) - A) %*% (diag(n) - A) %*% Y) -
          2 * log(max(0, 1 - tr(A) / n - 2 / n))
      })
    }
    else if (mode == "gmpml"){
      CV <- sapply(lambda, function(k){
        A <- K_mat %*% ginv(K_mat + k * diag(n))
        log(t(Y) %*% (diag(n) - A) %*% Y) -
          1 / (n - 1) * log(det((diag(n) - A)))
      })
    }
    else
      stop("mode must be loocv, AICc, GCVc or gmpml!")

    lambda0 <- lambda[which(CV == min(CV))]
    return(lambda0)
  }
