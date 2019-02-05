#' Calculating Tuning Parameters
#' 
#' Calculate tuning parameters based on given criteria.
#' 
#' There are seven tuning parameter selections here:
#' 
#' \bold{leave-one-out Cross Validation}
#' 
#' \deqn{\lambda_{n-CV}=\underset{\lambda \in
#' \Lambda}{argmin}\;\Big\{log\;y^{\star
#' T}[I-diag(A_\lambda)-\frac{1}{n}I]^{-1}(I-A_\lambda)^2[I-diag(A_\lambda)-
#' \frac{1}{n}I]^{-1}y^\star \Big\}}
#' 
#' \bold{Akaike Information Criteria}
#' 
#' \deqn{\lambda_{AIC}=\underset{\lambda \in \Lambda}{argmin}\Big\{log\;
#' y^{\star T}(I-A_\lambda)^2y^\star+\frac{2[tr(A_\lambda)+2]}{n}\Big\}}
#' 
#' \bold{Akaike Information Criteria (small sample size)}
#' 
#' \deqn{\lambda_{AICc}=\underset{\lambda \in \Lambda}{argmin}\Big\{log\;
#' y^{\star
#' T}(I-A_\lambda)^2y^\star+\frac{2[tr(A_\lambda)+2]}{n-tr(A_\lambda)-3}\Big\}}
#' 
#' \bold{Bayesian Information Criteria}
#' 
#' \deqn{\lambda_{BIC}=\underset{\lambda \in \Lambda}{argmin}\Big\{log\;
#' y^{\star T}(I-A_\lambda)^2y^\star+\frac{log(n)[tr(A_\lambda)+2]}{n}\Big\}}
#' 
#' \bold{Generalized Cross Validation}
#' 
#' \deqn{\lambda_{GCV}=\underset{\lambda \in \Lambda}{argmin}\Big\{log\;
#' y^{\star
#' T}(I-A_\lambda)^2y^\star-2log[1-\frac{tr(A_\lambda)}{n}-\frac{1}{n}]_+\Big\}}
#' 
#' \bold{Generalized Cross Validation (small sample size)}
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
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X (dataframe, n*p) Fixed effects variables in the dataframe (could
#' contains several subfactors).
#' @param K_mat (matrix, n*n) Estimated ensemble kernel matrix.
#' @param mode (character) A character string indicating which tuning parameter
#' criteria is to be used.
#' @param lambda (numeric) A numeric string specifying the range of tuning parameter 
#' to be chosen. The lower limit of lambda must be above 0.
#' @return \item{lambda0}{(numeric) The estimated tuning parameter.}
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
#' 
#' lambda0 <- tuning(Y = CVEK:::model_matrices$y, 
#' X = CVEK:::model_matrices$X, K_mat = CVEK:::K_ens, 
#' mode = "loocv", lambda = exp(seq(-10, 5)))
#' 
#' 
#' 
#' @export tuning
tuning <-
  function(Y, X, K_mat, mode, lambda) {
    
    mode <- match.arg(mode, c("AIC", "AICc", "BIC", "GCV", "GCVc", "gmpml", "loocv"))
    func_name <- paste0("tuning_", mode)
    lambda_selected <- do.call(func_name, 
                               list(Y = Y, X = X, K_mat = K_mat, lambda = lambda))
    if (length(lambda_selected) != 1) {
      warning(paste0("Multiple (", length(lambda_selected), 
                     ") optimal lambda's found, returning the smallest one."))
    }
    
    min(lambda_selected)
  }




#' Calculating Tuning Parameters Using AIC
#' 
#' Calculate tuning parameters based on AIC.
#' 
#' \bold{Akaike Information Criteria}
#' 
#' \deqn{\lambda_{AIC}=\underset{\lambda \in \Lambda}{argmin}\Big\{log\;
#' y^{\star T}(I-A_\lambda)^2y^\star+\frac{2[tr(A_\lambda)+2]}{n}\Big\}}
#' 
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X (dataframe, n*p) Fixed effects variables in the dataframe (could
#' contains several subfactors).
#' @param K_mat (matrix, n*n) Estimated ensemble kernel matrix.
#' @param lambda (numeric) A numeric string specifying the range of tuning parameter 
#' to be chosen. The lower limit of lambda must be above 0.
#' @return \item{lambda0}{(numeric) The estimated tuning parameter.}
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
tuning_AIC <-
  function(Y, X, K_mat, lambda) {
    
    n <- length(Y)
    CV <- sapply(lambda, function(k) {
      
      proj_matrix <- 
        estimate_ridge(X = X, K = K_mat, Y = Y, lambda = k)$proj_matrix
      A <- proj_matrix$total
      log(t(Y) %*% (diag(n) - A) %*% (diag(n) - A) %*% Y) +
        2 * (tr(A) + 2) / n
    })
    
    lambda[which(CV == min(CV))]
  }




#' Calculating Tuning Parameters Using AICc
#' 
#' Calculate tuning parameters based on AICc.
#' 
#' \bold{Akaike Information Criteria (small sample size)}
#' 
#' \deqn{\lambda_{AICc}=\underset{\lambda \in \Lambda}{argmin}\Big\{log\;
#' y^{\star
#' T}(I-A_\lambda)^2y^\star+\frac{2[tr(A_\lambda)+2]}{n-tr(A_\lambda)-3}\Big\}}
#' 
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X (dataframe, n*p) Fixed effects variables in the dataframe (could
#' contains several subfactors).
#' @param K_mat (matrix, n*n) Estimated ensemble kernel matrix.
#' @param lambda (numeric) A numeric string specifying the range of tuning parameter 
#' to be chosen. The lower limit of lambda must be above 0.
#' @return \item{lambda0}{(numeric) The estimated tuning parameter.}
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
tuning_AICc <-
  function(Y, X, K_mat, lambda) {
    
    n <- length(Y)
    CV <- sapply(lambda, function(k) {

      proj_matrix <- 
        estimate_ridge(X = X, K = K_mat, Y = Y, lambda = k)$proj_matrix
      A <- proj_matrix$total
      log(t(Y) %*% (diag(n) - A) %*% (diag(n) - A) %*% Y) +
        2 * (tr(A) + 2) / (n - tr(A) - 3)
    })
    
    lambda[which(CV == min(CV))]
  }




#' Calculating Tuning Parameters Using BIC
#' 
#' Calculate tuning parameters based on BIC.
#' 
#' \bold{Bayesian Information Criteria}
#' 
#' \deqn{\lambda_{BIC}=\underset{\lambda \in \Lambda}{argmin}\Big\{log\;
#' y^{\star T}(I-A_\lambda)^2y^\star+\frac{log(n)[tr(A_\lambda)+2]}{n}\Big\}}
#' 
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X (dataframe, n*p) Fixed effects variables in the dataframe (could
#' contains several subfactors).
#' @param K_mat (matrix, n*n) Estimated ensemble kernel matrix.
#' @param lambda (numeric) A numeric string specifying the range of tuning parameter 
#' to be chosen. The lower limit of lambda must be above 0.
#' @return \item{lambda0}{(numeric) The estimated tuning parameter.}
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
tuning_BIC <-
  function(Y, X, K_mat, lambda) {
    
    n <- length(Y)
    CV <- sapply(lambda, function(k) {
      
      proj_matrix <- 
        estimate_ridge(X = X, K = K_mat, Y = Y, lambda = k)$proj_matrix
      A <- proj_matrix$total
      log(t(Y) %*% (diag(n) - A) %*% (diag(n) - A) %*% Y) +
        log(n) * (tr(A) + 2) / n
    })
    
    lambda[which(CV == min(CV))]
  }




#' Calculating Tuning Parameters Using GCV
#' 
#' Calculate tuning parameters based on GCV.
#' 
#' \bold{Generalized Cross Validation}
#' 
#' \deqn{\lambda_{GCV}=\underset{\lambda \in \Lambda}{argmin}\Big\{log\;
#' y^{\star
#' T}(I-A_\lambda)^2y^\star-2log[1-\frac{tr(A_\lambda)}{n}-\frac{1}{n}]_+\Big\}}
#' 
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X (dataframe, n*p) Fixed effects variables in the dataframe (could
#' contains several subfactors).
#' @param K_mat (matrix, n*n) Estimated ensemble kernel matrix.
#' @param lambda (numeric) A numeric string specifying the range of tuning parameter 
#' to be chosen. The lower limit of lambda must be above 0.
#' @return \item{lambda0}{(numeric) The estimated tuning parameter.}
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
tuning_GCV <-
  function(Y, X, K_mat, lambda) {
    
    n <- length(Y)
    CV <- sapply(lambda, function(k) {
      
      proj_matrix <- 
        estimate_ridge(X = X, K = K_mat, Y = Y, lambda = k)$proj_matrix      
      A <- proj_matrix$total
      
      log(t(Y) %*% (diag(n) - A) %*% (diag(n) - A) %*% Y) -
        2 * log(1 - tr(A) / n - 1 / n)
    })
    
    lambda[which(CV == min(CV))]
  }





#' Calculating Tuning Parameters Using GCVc
#' 
#' Calculate tuning parameters based on GCVc.
#' 
#' \bold{Generalized Cross Validation (small sample size)}
#' 
#' \deqn{\lambda_{GCVc}=\underset{\lambda \in \Lambda}{argmin}\Big\{log\;
#' y^{\star
#' T}(I-A_\lambda)^2y^\star-2log[1-\frac{tr(A_\lambda)}{n}-\frac{2}{n}]_+\Big\}}
#' 
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X (dataframe, n*p) Fixed effects variables in the dataframe (could
#' contains several subfactors).
#' @param K_mat (matrix, n*n) Estimated ensemble kernel matrix.
#' @param lambda (numeric) A numeric string specifying the range of tuning parameter 
#' to be chosen. The lower limit of lambda must be above 0.
#' @return \item{lambda0}{(numeric) The estimated tuning parameter.}
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
tuning_GCVc <-
  function(Y, X, K_mat, lambda) {
    
    n <- length(Y)
    CV <- sapply(lambda, function(k) {
      
      proj_matrix <- 
        estimate_ridge(X = X, K = K_mat, Y = Y, lambda = k)$proj_matrix      
      A <- proj_matrix$total
      
      log(t(Y) %*% (diag(n) - A) %*% (diag(n) - A) %*% Y) -
        2 * log(max(0, 1 - tr(A) / n - 2 / n))
    })
    
    lambda[which(CV == min(CV))]
  }





#' Calculating Tuning Parameters Using GMPML
#' 
#' Calculate tuning parameters based on Generalized Maximum Profile Marginal
#' Likelihood.
#' 
#' \bold{Generalized Maximum Profile Marginal Likelihood}
#' 
#' \deqn{\lambda_{GMPML}=\underset{\lambda \in \Lambda}{argmin}\Big\{log\;
#' y^{\star T}(I-A_\lambda)y^\star-\frac{1}{n-1}log \mid I-A_\lambda \mid
#' \Big\}}
#' 
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X (dataframe, n*p) Fixed effects variables in the dataframe (could
#' contains several subfactors).
#' @param K_mat (matrix, n*n) Estimated ensemble kernel matrix.
#' @param lambda (numeric) A numeric string specifying the range of tuning parameter 
#' to be chosen. The lower limit of lambda must be above 0.
#' @return \item{lambda0}{(numeric) The estimated tuning parameter.}
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
tuning_gmpml <-
  function(Y, X, K_mat, lambda) {
    
    n <- length(Y)
    CV <- sapply(lambda, function(k){
      
      A_kernel_only <- K_mat %*% ginv(K_mat + k * diag(n))
      log_det <- unlist(determinant(diag(n) - A_kernel_only), 
                        use.names = FALSE)[1]
      proj_matrix <- 
        estimate_ridge(X = X, K = K_mat, Y = Y, lambda = k)$proj_matrix      
      A <- proj_matrix$total
      log(t(Y) %*% (diag(n) - A) %*% Y) - 1 / (n - 1) * log_det
    })
    
    lambda[which(CV == min(CV))]
  }





#' Calculating Tuning Parameters Using looCV
#' 
#' Calculate tuning parameters based on given leave-one-out Cross Validation.
#' 
#' \bold{leave-one-out Cross Validation}
#' 
#' \deqn{\lambda_{n-CV}=\underset{\lambda \in
#' \Lambda}{argmin}\;\Big\{log\;y^{\star
#' T}[I-diag(A_\lambda)-\frac{1}{n}I]^{-1}(I-A_\lambda)^2[I-diag(A_\lambda)-
#' \frac{1}{n}I]^{-1}y^\star \Big\}}
#' 
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X (dataframe, n*p) Fixed effects variables in the dataframe (could
#' contains several subfactors).
#' @param K_mat (matrix, n*n) Estimated ensemble kernel matrix.
#' @param lambda (numeric) A numeric string specifying the range of tuning parameter 
#' to be chosen. The lower limit of lambda must be above 0.
#' @return \item{lambda0}{(numeric) The estimated tuning parameter.}
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
tuning_loocv <-
  function(Y, X, K_mat, lambda) {
    
    n <- length(Y)
    CV <- sapply(lambda, function(k) {
      proj_matrix <- 
        estimate_ridge(X = X, K = K_mat, Y = Y, lambda = k)$proj_matrix      
      A <- proj_matrix$total
      
      sum(((diag(n) - A) %*% Y / diag(diag(n) - A)) ^ 2)
    })
    
    lambda[which(CV == min(CV))]
  }

