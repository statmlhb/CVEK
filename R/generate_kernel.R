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
#' x')=\frac{2}{\pi}sin^{-1}\Big(\frac{2\tilde{x}^T
#' \tilde{x}'}{\sqrt{(1+2\tilde{x}^T \tilde{x})(1+2\tilde{x}'^T
#' \tilde{x}')}}\Big)}
#'
#' @param method (character) A character string indicating which kernel 
#' is to be computed.
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
#'                                     kern_par[d, ]$l,
#'                                     kern_par[d, ]$p)
#' }
#'
#'
#' @export generate_kernel
#'
#' @import mvtnorm MASS psych limSolve stats

generate_kernel <-
  function(method = "rbf", l = 1, p = 2) {
    method <- match.arg(method, c("intercept", "linear", 
                                  "polynomial", "rbf", 
                                  "matern", "rational", "nn"))
    func_name <- paste0("kernel_", method)
    vector_wise <- do.call(func_name, list(l = l, p = p))
    
    kern <- function(X2, X1) 
      apply(X1, 1, function(xp){
        apply(X2, 1, function(xq){
          vector_wise(xp, xq, l, p)
        }
        )
      }
      )
    
    kern
  }


#' Generating A Single Vector-wise Function Using Intercept
#'
#' Generate vector-wise functions for two vectors using intercept kernel.
#'
#' \bold{Polynomial Kernels} \deqn{k(x, x')=(x \cdot x')^p} We have intercept
#' kernel when \eqn{p=0}, and linear kernel when \eqn{p=1}.
#'
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{vector_wise}{(function) A function calculating 
#' the relevance of two vectors.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#'
#' @export kernel_intercept
kernel_intercept <-
  function(l, p) {
    vector_wise <- function(xp, xq, l, p) {
      1
    }
    vector_wise
  }


#' Generating A Single Vector-wise Function Using Linear
#'
#' Generate vector-wise functions for two vectors using linear kernel.
#'
#' \bold{Polynomial Kernels} \deqn{k(x, x')=(x \cdot x')^p} We have intercept
#' kernel when \eqn{p=0}, and linear kernel when \eqn{p=1}.
#'
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{vector_wise}{(function) A function calculating 
#' the relevance of two vectors.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#'
#' @export kernel_linear
kernel_linear <-
  function(l, p) {
    vector_wise <- function(xp, xq, l, p) {
      crossprod(xp, xq)
    }
    vector_wise
  }


#' Generating A Single Vector-wise Function Using Polynomial
#'
#' Generate vector-wise functions for two vectors using polynomial kernel.
#'
#' \bold{Polynomial Kernels} \deqn{k(x, x')=(x \cdot x')^p} We have intercept
#' kernel when \eqn{p=0}, and linear kernel when \eqn{p=1}.
#'
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{vector_wise}{(function) A function calculating 
#' the relevance of two vectors.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#'
#' @export kernel_polynomial
kernel_polynomial <-
  function(l, p) {
    vector_wise <- function(xp, xq, l, p) {
      (crossprod(xp, xq) + 1) ^ p
    }
    vector_wise
  }


#' Generating A Single Vector-wise Function Using RBF
#'
#' Generate vector-wise functions for two vectors using rbf kernel.
#'
#' \bold{Gaussian RBF Kernels} \deqn{k_{SE}(r)=exp\Big(-\frac{r^2}{2l^2}\Big)}
#'
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{vector_wise}{(function) A function calculating 
#' the relevance of two vectors.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#'
#' @export kernel_rbf
kernel_rbf <-
  function(l, p) {
    vector_wise <- function(xp, xq, l, p) {
      exp((2 * crossprod(xp, xq) - crossprod(xp) - crossprod(xq)) / (2 * l ^ 2))
    }
    vector_wise
  }


#' Generating A Single Vector-wise Function Using Matern
#'
#' Generate vector-wise functions for two vectors using matern kernel.
#'
#' \bold{Matern Kernels}
#' \deqn{k_{Matern}(r)=\frac{2^{1-\nu}}{\Gamma(\nu)}\Big(\frac{\sqrt{2\nu
#' r}}{l}\Big)^\nu K_\nu \Big(\frac{\sqrt{2\nu r}}{l}\Big)}
#'
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{vector_wise}{(function) A function calculating 
#' the relevance of two vectors.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#'
#' @export kernel_matern
kernel_matern <-
  function(l, p) {
    vector_wise <- function(xp, xq, l, p){
      r <- sqrt(crossprod(xp) + crossprod(xq) - 2 * crossprod(xp, xq))
      v <- p + 1 / 2
      s <- 0
      for (i in 0:p) {
        s <- s + factorial(p + i) / (factorial(i) * factorial(p - i)) *
          (sqrt(8 * v) * r / l) ^ (p - i)
      }
      exp(-sqrt(2 * v) * r / l) * gamma(p + 1) / gamma(2 * p + 1) * s
    }
    vector_wise
  }


#' Generating A Single Vector-wise Function Using Rational Quadratic
#'
#' Generate vector-wise functions for two vectors using rational kernel.
#'
#' \bold{Rational Quadratic Kernels} \deqn{k_{RQ}(r)=\Big(1+\frac{r^2}{2\alpha
#' l^2}\Big)^{-\alpha}}
#'
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{vector_wise}{(function) A function calculating 
#' the relevance of two vectors.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#'
#' @export kernel_rational
kernel_rational <-
  function(l, p) {
    vector_wise <- function(xp, xq, l, p){
      r <- sqrt(crossprod(xp) + crossprod(xq) - 2 * crossprod(xp, xq))
      (1 + r ^ 2 / (2 * p * l ^ 2)) ^ (- p)
    }
    vector_wise
  }


#' Generating A Single Vector-wise Function Using Neural Network
#'
#' Generate vector-wise functions for two vectors using neural network kernel.
#'
#' \bold{Neural Network Kernels} \deqn{k_{NN}(x,
#' x')=\frac{2}{\pi}sin^{-1}\Big(\frac{2\tilde{x}^T
#' \tilde{x}'}{\sqrt{(1+2\tilde{x}^T \tilde{x})(1+2\tilde{x}'^T
#' \tilde{x}')}}\Big)}
#'
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{vector_wise}{(function) A function calculating 
#' the relevance of two vectors.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#'
#' @export kernel_nn
kernel_nn <-
  function(l, p) {
    vector_wise <- function(xp, xq, l, p){
      xp <- c(1, xp)
      xq <- c(1, xq)
      s <- 2 * crossprod(xp, xq) / (sqrt((1 + 2 * crossprod(xp, xp))
                                         * (1 + 2 * crossprod(xq, xq))))
      asin(s)
    }
    vector_wise
  }

