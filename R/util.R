#' From Vectors to Single Variables
#' 
#' Transform format of predictors from vectors to single variables.
#' 
#' 
#' @param formula (formula) A symbolic description of the model to be fitted.
#' @param label_names (list) A character string indicating all the interior
#' variables included in each group of random effect.
#' @return \item{generic_formula}{(formula) A symbolic description of the model
#' written in single variables format.}
#' 
#' \item{length_main}{(integer) A numeric value indicating the length of main
#' random effects.}
#' @author Wenying Deng
#' @examples
#' 
#' 
#' 
#' generic_formula0 <- generate_formula(formula = Y ~ X + Z1 + Z2,
#' label_names = list(Z1 = c("z1", "z2"), Z2 = c("z3", "z4")))
#' 
#' 
#' 
#' @export generate_formula
generate_formula <-
  function(formula, label_names) {
    
    formula_factors <- attr(terms(formula), "factors")
    generic_formula <- Y ~ 1
    length_main <- 0
    for (i in 2:dim(formula_factors)[2]) {
      terms_names <-
        rownames(formula_factors)[which(formula_factors[, i] == 1)]
      if (length(terms_names) == 1) {
        generic_formula <-
          update.formula(generic_formula,
                         as.formula(paste(as.character(
                           attr(terms(formula), "variables"))[2],
                           paste(label_names[[terms_names]], 
                                 collapse=" + "), sep=" ~ .+")))
        length_main <- length_main + length(label_names[[terms_names]])
      } else {
        interaction_formula <-
          paste("(", paste(label_names[[terms_names[1]]], 
                           collapse=" + "), ")", sep="")
        for (j in 2:length(terms_names)) {
          interaction_formula <-
            paste(interaction_formula, "*(",
                  paste(label_names[[terms_names[j]]], 
                        collapse=" + "), ")", sep=" ")
        }
        generic_formula <-
          update.formula(generic_formula,
                         as.formula(paste(as.character(
                           attr(terms(formula), "variables"))[2], 
                           interaction_formula, sep=" ~ .+")
                         )
          )
      }
    }
    
    list(generic_formula =  generic_formula, length_main = length_main)
  }






#' Generating Original Data
#' 
#' Generate original data based on specific kernels.
#' 
#' This function generates with a specific dataset. The argument int_effect
#' represents the strength of interaction of random effects relative to the
#' main random effects since all sampled functions have been standardized to
#' have unit norm.
#' 
#' @param n (integer) A numeric number specifying the number of observations.
#' @param fixed_num (integer) A numeric number specifying the dimension of
#' fixed effects.
#' @param label_names (list) A character string indicating all the interior
#' variables included in each group of random effects.
#' @param method (character) A character string indicating which kernel is to
#' be computed.
#' @param l (numeric) A numeric number indicating the hyperparameter
#' (flexibility) of a specific kernel.
#' @param d (integer) For polynomial, d is the power; for matern, v = d + 1 /
#' 2; for rational, alpha = d.
#' @param int_effect (numeric) A numeric number specifying the size of
#' interaction.
#' @param eps (numeric) A numeric number indicating the size of noise of fixed
#' effects.
#' @return \item{data}{(dataframe, n*(p+q)) A dataframe to be fitted.}
#' @author Wenying Deng
#' @examples
#' 
#' 
#' 
#' mydata <- generate_data(n = 100, fixed_num = 1, label_names =
#' list(Z1 = c("z1", "z2"), Z2 = c("z3", "z4")),
#' method = "rbf", l = 1, d = 2, int_effect = 0, eps = .01)
#' 
#' 
#' 
#' @export generate_data
generate_data <-
  function(n, fixed_num = 1, 
           label_names = NULL, 
           method = "rbf", 
           l = 1, d = 2, 
           int_effect = 0, eps = .01) {
    
    if ((fixed_num == 0) & is.null(label_names)) {
      stop("fixed effect and random effect can not be null simultaneously!")
    }
    if (!is.null(label_names)) {
      Z1 <- rmvnorm(n = n,
                    mean = rep(0, length(label_names[[1]])),
                    sigma = diag(length(label_names[[1]])))
      Z2 <- rmvnorm(n = n,
                    mean = rep(0, length(label_names[[2]])),
                    sigma = diag(length(label_names[[2]])))
      kern <- generate_kernel(method = method, l = l, d = d)
      w1 <- rnorm(n)
      w2 <- w1
      w12 <- rnorm(n)
      K1 <- kern(Z1, Z1)
      K2 <- kern(Z2, Z2)
      K1 <- K1 / tr(K1)
      K2 <- K2 / tr(K2)
      h0 <- K1 %*% w1 + K2 %*% w2
      h0 <- h0 / sqrt(sum(h0 ^ 2))
      h1_prime <- (K1 * K2) %*% w12
      Ks <- svd(K1 + K2)
      if (length(Ks$d / sum(Ks$d) > .001) > 0) {
        len <- length(Ks$d[Ks$d / sum(Ks$d) > .001])
        U0 <- Ks$u[, 1:len]
        h1_prime_hat <- fitted(lm(h1_prime ~ U0))
        h1 <- h1_prime - h1_prime_hat
        if (all(h1 == 0)) {
          warning("interaction term colinear with main-effect space!")
          h1 <- h1_prime
          h1 <- h1 / sqrt(sum(h1 ^ 2))
        } else {
          h1 <- h1 / sqrt(sum(h1 ^ 2))
        }
      } else {
        warning("largest eigen value smaller than 0.001!")
        h1 <- h1_prime
        h1 <- h1 / sqrt(sum(h1 ^ 2))
      }
    } else {
      h0 <- 0
      h1 <- 0
    }
    if (fixed_num > 0) {
      X <- rmvnorm(n = n,
                   mean = rep(0, fixed_num),
                   sigma = diag(fixed_num))
      beta <- rnorm(fixed_num)
      Y_fixed <- X %*% beta
      Xnam <- paste0("x", 1:fixed_num)
    } else {
      X <- NULL
      Y_fixed <- 0
      Xnam <- NULL
    }
    
    Y <- Y_fixed + h0 + int_effect * h1 + rnorm(1) + rnorm(n, 0, eps)
    data <- as.data.frame(cbind(Y, X, Z1, Z2))
    colnames(data) <- c("Y", Xnam, label_names[[1]], label_names[[2]])
    
    data
  }





#' Estimating Noise
#' 
#' An implementation of Gaussian processes for estimating noise.
#' 
#' 
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X (dataframe, n*p) Fixed effects variables in the dataframe (could
#' contains several subfactors).
#' @param lambda_hat (numeric) The selected tuning parameter based on the
#' estimated ensemble kernel matrix.
#' @param y_fixed_hat (vector of length n) Estimated fixed effects of the
#' responses.
#' @param alpha_hat (vector of length n) Random effects estimator of the
#' estimated ensemble kernel matrix.
#' @param K_hat (matrix, n*n) Estimated ensemble kernel matrix.
#' @return \item{sigma2_hat}{(numeric) The estimated noise of the fixed
#' effects.}
#' @author Wenying Deng
#' @references Jeremiah Zhe Liu and Brent Coull. Robust Hypothesis Test for
#' Nonlinear Effect with Gaussian Processes. October 2017.s
#' @keywords internal
#' @export estimate_sigma2
estimate_sigma2 <- function(Y, X, lambda_hat, y_fixed_hat, alpha_hat, K_hat) {

  n <- length(Y)
  V_inv <- ginv(K_hat + lambda_hat * diag(n))
  B_mat <- ginv(t(X) %*% V_inv %*% X) %*% t(X) %*% V_inv
  P_X <- X %*% B_mat
  P_K <- K_hat %*% V_inv %*% (diag(n) - P_X)
  A <- P_X + P_K
  sigma2_hat <- sum((Y - y_fixed_hat - K_hat %*% alpha_hat) ^ 2) / (n - tr(A) - 1)
  
  sigma2_hat
}



#' Computing Score Test Statistics.
#' 
#' Compute score test statistics.
#' 
#' The test statistic is distributed as a scaled Chi-squared distribution.
#' 
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param K_int (matrix, n*n) The kernel matrix to be tested.
#' @param y_fixed (vector of length n) Estimated fixed effects of the
#' responses.
#' @param K_0 (matrix, n*n) Estimated ensemble kernel matrix.
#' @param sigma2_hat (numeric) The estimated noise of the fixed effects.
#' @param tau_hat (numeric) The estimated noise of the random effects.
#' @return \item{test_stat}{(numeric) The computed test statistic.}
#' @author Wenying Deng
#' @references Arnab Maity and Xihong Lin. Powerful tests for detecting a gene
#' effect in the presence of possible gene-gene interactions using garrote
#' kernel machines. December 2011.
#' @keywords internal
#' @export compute_stat
compute_stat <-
  function(Y, K_int, y_fixed, K0, sigma2_hat, tau_hat) {
    
    n <- length(Y)
    V0_inv <- ginv(tau_hat * K0 + sigma2_hat * diag(n))
    test_stat <- tau_hat * t(Y - y_fixed) %*% V0_inv %*%
      K_int %*% V0_inv %*% (Y - y_fixed) / 2

    test_stat
  }



#' Computing Information Matrices
#' 
#' Compute information matrices based on block matrices.
#' 
#' This function gives the information value of the interaction strength.
#' 
#' @param P0_mat (matrix, n*n) Scale projection matrix under REML.
#' @param mat_del (matrix, n*n) Derivative of the scale covariance matrix of Y
#' with respect to delta.
#' @param mat_sigma2 (matrix, n*n) Derivative of the scale covariance matrix of
#' Y with respect to sigma2.
#' @param mat_tau (matrix, n*n) Derivative of the scale covariance matrix of Y
#' with respect to tau.
#' @return \item{I0}{(matrix, n*n) The computed information value.}
#' @author Wenying Deng
#' @references Arnab Maity and Xihong Lin. Powerful tests for detecting a gene
#' effect in the presence of possible gene-gene interactions using garrote
#' kernel machines. December 2011.
#' @keywords internal
#' @export compute_info
compute_info <-
  function(P0_mat, mat_del = NULL, mat_sigma2 = NULL, mat_tau = NULL) {
    
    I0 <- matrix(NA, 3, 3)
    I0[1, 1] <- tr(P0_mat %*% mat_del %*% P0_mat %*% mat_del) / 2  
    I0[1, 2] <- tr(P0_mat %*% mat_del %*% P0_mat %*% mat_sigma2) / 2
    I0[2, 1] <- I0[1, 2]
    I0[1, 3] <- tr(P0_mat %*% mat_del %*% P0_mat %*% mat_tau) / 2
    I0[3, 1] <- I0[1, 3]
    I0[2, 2] <- tr(P0_mat %*% mat_sigma2 %*% P0_mat %*% mat_sigma2) / 2
    I0[2, 3] <-  tr(P0_mat %*% mat_sigma2 %*% P0_mat %*% mat_tau) / 2  
    I0[3, 2] <- I0[2, 3]
    I0[3, 3] <-  tr(P0_mat %*% mat_tau %*% P0_mat %*% mat_tau) / 2  

    I0
  }



#' Standardizing Matrix
#' 
#' Center and scale the data matrix into mean zero and standard deviation one.
#' 
#' This function gives the standardized data matrix.
#' 
#' @param X (matrix, n*p0) Original data matrix.
#' @return \item{X}{(matrix, n*p0) Standardized data matrix.}
#' @author Wenying Deng
#' @keywords internal
#' @export standardize
standardize <- function(X) {
  
  Xm <- colMeans(X)
  n <- nrow(X)
  p <- ncol(X)
  X <- X - rep(Xm, rep(n, p))
  Xscale <- drop(rep(1 / n, n) %*% X ^ 2) ^ .5
  X <- X / rep(Xscale, rep(n, p))
  
  X
}



#' Computing Euclidean Distance between Two Vectors (Matrices)
#' 
#' Compute the L2 distance between two vectors or matrices.
#' 
#' This function gives the Euclidean distance between two 
#' vectors or matrices.
#' 
#' @param x1 (vector/matrix) The first vector/matrix.
#' @param x2 (vector/matrix, the same dimension as x1) 
#' The second vector/matrix.
#' @return \item{dist}{(numeric) Euclidean distance.}
#' @author Wenying Deng
#' @keywords internal
#' @export euc_dist
euc_dist <- function(x1, x2 = NULL) {
  if (is.null(x2)) {
    dist <- sqrt(sum(x1 ^ 2))
  } else {
    dist <- sqrt(sum((x1 - x2) ^ 2))
  }
  dist
}
