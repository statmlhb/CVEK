#' Predicting New Response
#' 
#' Predicting new response based on given design matrix and 
#' the estimation result.
#' 
#' After we obtain the estimation result, we can predict new response.
#' 
#' @param object (list) Estimation results returned by estimation() procedure.
#' @param data_new (dataframe) The new set of predictors, whose name is 
#' the same as those of formula estimation().
#' @return \item{y_pred}{(vector of length n) Predicted new response.}
#' @author Wenying Deng
#' @export predict.cvek
predict.cvek <- function(object, data_new) {
  
  model_matrices <- object$model_matrices
  kern_func_list <- object$kern_func_list
  new_matrices <- parse_cvek_formula(object$formula, 
                                     kern_func_list, 
                                     data = object$data, 
                                     data_new = data_new)
  X <- new_matrices$X
  n <- nrow(data)
  A <- 0
  Xmat <- ginv(t(model_matrices$X) 
               %*% model_matrices$X) %*% t(model_matrices$X)
  H <- model_matrices$X %*% Xmat
  H_star <- X %*% Xmat
  n <- length(object$alpha)
  P_K_star <- list()
  P_X_star <- list()
  y_pred <- 0
  for (k in seq(length(kern_func_list))) {
    
    B_temp <- 0
    for (d in seq(length(new_matrices$K[[k]]))) {
      S_d_star <- new_matrices$K[[k]][[d]] %*% ginv(model_matrices$K[[k]][[d]] 
                                                    + object$base_est$lambda_list[[k]] * diag(n))
      B_temp <- B_temp + S_d_star %*% (diag(n) + object$base_est$A_proc_list[[k]][[d]])
    }
    B_star <- B_temp %*% (diag(n) - object$base_est$P_K_hat[[k]])
    P_K <- ginv(diag(n) - object$base_est$P_K_hat[[k]] %*% H) %*% 
      object$base_est$P_K_hat[[k]] %*% (diag(n) - H)
    P_K_star[[k]] <- B_star %*% (diag(n) - H + H %*% object$base_est$P_K_hat[[k]])
    P_X_star[[k]] <- H_star %*% (diag(n) - P_K)
    y_pred <- y_pred + object$u_hat[k] * (P_K_star[[k]] + P_X_star[[k]]) %*% object$data$y
  }

  y_pred
}
