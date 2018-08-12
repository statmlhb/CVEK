#' Generating A Single Kernel
#'
#' Generate kernels for the kernel library.
#'
#' There are seven kinds of kernel available here. For convenience, we define
#' \eqn{r=\mid x-x'\mid}.
#'
#' \bold{Gaussian RBF Kernels} \deqn{k_{SE}(r)=exp\Big(-\frac{r^2}{2l^2}\Big)}
#'
#' \bold{Matern Kernels}
#' \deqn{k_{Matern}(r)=\frac{2^{1-\nu}}{\Gamma(\nu)}\Big(\frac{\sqrt{2\nu
#' r}}{l}\Big)^\nu K_\nu \Big(\frac{\sqrt{2\nu r}}{l}\Big)}
#'
#' \bold{Rational Quadratic Kernels} \deqn{k_{RQ}(r)=\Big(1+\frac{r^2}{2\alpha
#' l^2}\Big)^{-\alpha}}
#'
#' \bold{Polynomial Kernels} \deqn{k(x, x')=(x \cdot x')^p} We have intercept
#' kernel when \eqn{p=0}, and linear kernel when \eqn{p=1}.
#'
#' \bold{Neural Network Kernels} \deqn{k_{NN}(x,
#' x')=\frac{2}{\pi}sin^{-1}\Big(\frac{2\tilde{x}^T\Sigma
#' \tilde{x}'}{\sqrt{(1+2\tilde{x}^T\Sigma \tilde{x})(1+2\tilde{x}'^T\Sigma
#' \tilde{x}')}}\Big)}
#'
#' @param method (character) A character string indicating which kernel 
#' is to be computed.
#' @param Sigma (matrix) The covariance matrix for neural network kernel.
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{kern}{(function) A function indicating the generated kernel.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#' @examples
#'
#'
#' kern_list <- list()
#' for (d in 1:nrow(kern_par)) {
#'   kern_list[[d]] <- generate_kernel(kern_par[d, ]$method,
#'                                     kern_par[d, ]$Sigma,
#'                                     kern_par[d, ]$l,
#'                                     kern_par[d, ]$p)
#' }
#'
#'
#' @export generate_kernel
#'
#' @import mvtnorm MASS psych limSolve stats

generate_kernel <-
  function(method = "rbf", Sigma = 0, l = 1, p = 2) {
    method <- match.arg(method, c("intercept", "linear", 
                                  "polynomial", "rbf", 
                                  "matern", "rational", "nn"))
    func_name <- paste0("kernel_", method)
    point_wise <- do.call(func_name, list(Sigma = Sigma, l = l, p = p))
    
    kern <- function(X2, X1) 
      apply(X1, 1, function(xp){
        apply(X2, 1, function(xq){
          point_wise(xp, xq, Sigma, l, p)
        }
        )
      }
      )
    
    kern
  }


#' Generating A Single Point-wise Function Using Intercept
#'
#' Generate point-wise functions for two vectors using intercept kernel.
#'
#' \bold{Polynomial Kernels} \deqn{k(x, x')=(x \cdot x')^p} We have intercept
#' kernel when \eqn{p=0}, and linear kernel when \eqn{p=1}.
#'
#' @param Sigma (matrix) The covariance matrix for neural network kernel.
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{point_wise}{(function) A function calculating 
#' the relevance of two vectors.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#'
#' @export kernel_intercept
kernel_intercept <-
  function(Sigma, l, p) {
    point_wise <- function(xp, xq, Sigma, l, p) {
      1
    }
    point_wise
  }


#' Generating A Single Point-wise Function Using Linear
#'
#' Generate point-wise functions for two vectors using linear kernel.
#'
#' \bold{Polynomial Kernels} \deqn{k(x, x')=(x \cdot x')^p} We have intercept
#' kernel when \eqn{p=0}, and linear kernel when \eqn{p=1}.
#'
#' @param Sigma (matrix) The covariance matrix for neural network kernel.
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{point_wise}{(function) A function calculating 
#' the relevance of two vectors.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#'
#' @export kernel_linear
kernel_linear <-
  function(Sigma, l, p) {
    point_wise <- function(xp, xq, Sigma, l, p) {
      t(xp) %*% xq
    }
    point_wise
  }


#' Generating A Single Point-wise Function Using Polynomial
#'
#' Generate point-wise functions for two vectors using polynomial kernel.
#'
#' \bold{Polynomial Kernels} \deqn{k(x, x')=(x \cdot x')^p} We have intercept
#' kernel when \eqn{p=0}, and linear kernel when \eqn{p=1}.
#'
#' @param Sigma (matrix) The covariance matrix for neural network kernel.
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{point_wise}{(function) A function calculating 
#' the relevance of two vectors.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#'
#' @export kernel_polynomial
kernel_polynomial <-
  function(Sigma, l, p) {
    point_wise <- function(xp, xq, Sigma, l, p) {
      (t(xp) %*% xq + 1) ^ p
    }
    point_wise
  }


#' Generating A Single Point-wise Function Using RBF
#'
#' Generate point-wise functions for two vectors using rbf kernel.
#'
#' \bold{Gaussian RBF Kernels} \deqn{k_{SE}(r)=exp\Big(-\frac{r^2}{2l^2}\Big)}
#'
#' @param Sigma (matrix) The covariance matrix for neural network kernel.
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{point_wise}{(function) A function calculating 
#' the relevance of two vectors.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#'
#' @export kernel_rbf
kernel_rbf <-
  function(Sigma, l, p) {
    point_wise <- function(xp, xq, Sigma, l, p) {
      exp(- sum((xp - xq) ^ 2) / (2 * l ^ 2))
    }
    point_wise
  }


#' Generating A Single Point-wise Function Using Matern
#'
#' Generate point-wise functions for two vectors using matern kernel.
#'
#' \bold{Matern Kernels}
#' \deqn{k_{Matern}(r)=\frac{2^{1-\nu}}{\Gamma(\nu)}\Big(\frac{\sqrt{2\nu
#' r}}{l}\Big)^\nu K_\nu \Big(\frac{\sqrt{2\nu r}}{l}\Big)}
#'
#' @param Sigma (matrix) The covariance matrix for neural network kernel.
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{point_wise}{(function) A function calculating 
#' the relevance of two vectors.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#'
#' @export kernel_matern
kernel_matern <-
  function(Sigma, l, p) {
    point_wise <- function(xp, xq, Sigma, l, p){
      r <- sqrt(sum((xp - xq) ^ 2))
      v <- p + 1 / 2
      s <- 0
      for (i in 0:p) {
        s <- s + factorial(p + i) / (factorial(i) * factorial(p - i)) *
          (sqrt(8 * v) * r / l) ^ (p - i)
      }
      exp(-sqrt(2 * v) * r / l) * gamma(p + 1) / gamma(2 * p + 1) * s
    }
    point_wise
  }


#' Generating A Single Point-wise Function Using Rational Quadratic
#'
#' Generate point-wise functions for two vectors using rational kernel.
#'
#' \bold{Rational Quadratic Kernels} \deqn{k_{RQ}(r)=\Big(1+\frac{r^2}{2\alpha
#' l^2}\Big)^{-\alpha}}
#'
#' @param Sigma (matrix) The covariance matrix for neural network kernel.
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{point_wise}{(function) A function calculating 
#' the relevance of two vectors.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#'
#' @export kernel_rational
kernel_rational <-
  function(Sigma, l, p) {
    point_wise <- function(xp, xq, Sigma, l, p){
      r <- sqrt(sum((xp - xq) ^ 2))
      (1 + r ^ 2 / (2 * p * l ^ 2)) ^ (- p)
    }
    point_wise
  }


#' Generating A Single Point-wise Function Using Neural Network
#'
#' Generate point-wise functions for two vectors using neural network kernel.
#'
#' \bold{Neural Network Kernels} \deqn{k_{NN}(x,
#' x')=\frac{2}{\pi}sin^{-1}\Big(\frac{2\tilde{x}^T\Sigma
#' \tilde{x}'}{\sqrt{(1+2\tilde{x}^T\Sigma \tilde{x})(1+2\tilde{x}'^T\Sigma
#' \tilde{x}')}}\Big)}
#'
#' @param Sigma (matrix) The covariance matrix for neural network kernel.
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{point_wise}{(function) A function calculating 
#' the relevance of two vectors.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#'
#' @export kernel_nn
kernel_nn <-
  function(Sigma, l, p) {
    point_wise <- function(xp, xq, Sigma, l, p){
      xp <- c(1, xp)
      xq <- c(1, xq)
      if(Sigma == 0) Sigma <- diag(length(xp))
      s <- 2 * t(xp) %*% Sigma %*% xq / (sqrt((1 + 2 * t(xp) %*% Sigma %*% xp)
                                              * (1 + 2 * t(xq) %*% Sigma %*% xq)))
      asin(s)
    }
    point_wise
  }

