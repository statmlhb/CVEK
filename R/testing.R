#' Conducting Score Tests for Interaction
#'
#' Conduct score tests comparing a fitted model and a more general alternative
#' model.
#'
#' There are two tests available here:
#'
#' \bold{Asymptotic Test}
#'
#' This is based on the classical variance component test to construct a
#' testing procedure for the hypothesis about Gaussian process function.
#'
#' \bold{Bootstrap Test}
#'
#' When it comes to small sample size, we can use bootstrap test instead, which
#' can give valid tests with moderate sample sizes and requires similar
#' computational effort to a permutation test.
#'
#' @param formula_int A symbolic description of the model with interaction.
#' @param label_names A character string indicating all the interior variables
#' included in each predictor.
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
#' @param test A character string indicating which test is to be used.
#' @param lambda A numeric string specifying the range of noise to be chosen.
#' The lower limit of lambda must be above 0.
#' @param B A numeric value indicating times of resampling when test = "boot".
#' @return \item{pvalue}{p-value of the test.}
#' @author Wenying Deng
#' @seealso method: \code{\link{kernelGenerate}}
#'
#' mode: \code{\link{tuning}}
#'
#' strategy: \code{\link{ensemble}}
#' @references Xihong Lin. Variance component testing in generalised linear
#' models with random effects. June 1997.
#'
#' Arnab Maity and Xihong Lin. Powerful tests for detecting a gene effect in
#' the presence of possible gene-gene interactions using garrote kernel
#' machines. December 2011.
#'
#' Petra Bu ̊zˇkova ́, Thomas Lumley, and Kenneth Rice. Permutation and
#' parametric bootstrap tests for gene-gene and gene-environment interactions.
#' January 2011.
#' @examples
#'
#'
#' ##testing(formula_int = Y ~ X1 * X2,
#' ##label_names = list(X1 = c("x1", "x2"), X2 = c("x3", "x4")),
#' ##Y, X1, X2, Kernlist, mode = "loocv", strategy = "erm",
#' ##beta = 1, test = "boot", lambda = exp(seq(-5, 5)), B = 100)
#'
#'
#' @export testing
testing <- function(formula_int, label_names, Y, X1, X2, Kernlist,
                    mode = "loocv", strategy = "erm", beta = 1,
                    test = "boot", lambda = exp(seq(-5, 5)), B = 100){

  re <- genericFormula(formula_int, label_names)
  generic_formula0 <- re$generic_formula
  len <- re$length_main
  data <- as.data.frame(cbind(Y, X1, X2))
  colnames(data) <- c("Y", label_names[[1]], label_names[[2]])
  X <- model.matrix(generic_formula0, data)[, -1]
  X12 <- X[, c((len + 1):dim(X)[2])]
  n <- length(Y)

  result <- estimation(Y, X1, X2, Kernlist, mode, strategy, beta, lambda)
  lam <- result[[1]]
  beta0 <- result[[2]]
  alpha0 <- result[[3]]
  K_gpr <- result[[4]]
  # u_weight <- result[[5]]
  sigma2_hat <- noiseEstimate(Y, lam, beta0, alpha0, K_gpr)
  tau_hat <- sigma2_hat / lam

  if (test == "boot"){
    # conduct bootstrap
    meanY <- K_gpr %*% alpha0 + beta0
    bs_test <- sapply(1:B, function(k){
      Ystar <- meanY + rnorm(n, sd = sqrt(sigma2_hat))
      scoreStat(n, Ystar, X12, beta0, sigma2_hat, tau_hat, K_gpr)
    })

    # assemble test statistic
    original_test <-
      scoreStat(n, Y, X12, beta0, sigma2_hat, tau_hat, K_gpr)

    pvalue <- sum(as.numeric(original_test) <= bs_test) / B
  }
  else if (test == "asym"){
    score_chi <-
      scoreStat(n, Y, X12, beta0, sigma2_hat, tau_hat, K_gpr)

    K0 <- K_gpr
    K12 <- X12 %*% t(X12)
    V0_inv <- ginv(tau_hat * K0 + sigma2_hat * diag(n))
    one <- rep(1, n)
    P0_mat <- V0_inv - V0_inv %*%
      one %*% ginv(t(one) %*% V0_inv %*% one) %*% t(one) %*% V0_inv

    drV0_tau <- K0
    drV0_sigma2 <- diag(n)
    drV0_del <- tau_hat * K12

    I0 <- infoMat(P0_mat,
                  mat_del = drV0_del, mat_sigma2 = drV0_sigma2,
                  mat_tau = drV0_tau)

    #Effective Info for delta
    tot_dim <- ncol(I0)
    I_deldel <-
      I0[1, 1] -
      I0[1, 2:tot_dim] %*% ginv(I0[2:tot_dim, 2:tot_dim]) %*% I0[2:tot_dim, 1]

    md <- tau_hat * tr(K12 %*% P0_mat) / 2

    m_chi <- I_deldel / (2 * md)
    d_chi <- md / m_chi

    pvalue <- 1 - pchisq(score_chi / m_chi, d_chi)
  }
  else
    stop("test must be boot or asym!")

  return(pvalue)
}
