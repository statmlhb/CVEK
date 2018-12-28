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
#' \bold{Polynomial Kernels} \deqn{k(x, x')=(x \cdot x')^d} We have intercept
#' kernel when \eqn{d=0}, and linear kernel when \eqn{d=1}.
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
#' @param d (integer) For polynomial, d is the power; for matern, v = d + 1 / 2; for
#' rational, alpha = d.
#' @return \item{kern}{(function) A function indicating the generated kernel.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#' @examples
#'
#'
#' kern_list <- list()
#' for (k in 1:nrow(kern_par)) {
#'   kern_list[[k]] <- generate_kernel(kern_par[k, ]$method,
#'                                     kern_par[k, ]$l,
#'                                     kern_par[k, ]$d)
#' }
#'
#'
#' @export generate_kernel
#'
#' @import mvtnorm MASS psych limSolve stats



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
#' \bold{Polynomial Kernels} \deqn{k(x, x')=(x \cdot x')^d} We have intercept
#' kernel when \eqn{d=0}, and linear kernel when \eqn{d=1}.
#' 
#' \bold{Neural Network Kernels} \deqn{k_{NN}(x,
#' x')=\frac{2}{\pi}sin^{-1}\Big(\frac{2\tilde{x}^T
#' \tilde{x}'}{\sqrt{(1+2\tilde{x}^T \tilde{x})(1+2\tilde{x}'^T
#' \tilde{x}')}}\Big)}
#' 
#' @param method (character) A character string indicating which kernel is to
#' be computed.
#' @param l (numeric) A numeric number indicating the hyperparameter
#' (flexibility) of a specific kernel.
#' @param d (integer) For polynomial, d is the power; for matern, v = d + 1 /
#' 2; for rational, alpha = d.
#' @return \item{kern}{(function) A function indicating the generated kernel.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#' @examples
#' 
#' 
#' 
#' kern_list <- list()
#' for (k in 1:nrow(kern_par)) {
#'   kern_list[[k]] <- generate_kernel(kern_par[k, ]$method,
#'                                     kern_par[k, ]$l,
#'                                     kern_par[k, ]$d)
#' }
#' 
#' 
#' 
#' @export generate_kernel
generate_kernel <-
  function(method = "rbf", l = 1, d = 2) {
    method <- match.arg(method, c("intercept", "linear", 
                                  "polynomial", "rbf", 
                                  "matern", "rational", "nn"))
    func_name <- paste0("kernel_", method)
    matrix_wise <- do.call(func_name, list(l = l, d = d))
    kern <- function(X1, X2) {
      matrix_wise(X1, X2, l, d)
    }
    
    kern
  }




#' Generating A Single Matrix-wise Function Using Intercept
#' 
#' Generate matrix-wise functions for two matrices using intercept kernel.
#' 
#' \bold{Polynomial Kernels} \deqn{k(x, x')=(x \cdot x')^d} We have intercept
#' kernel when \eqn{d=0}, and linear kernel when \eqn{d=1}.
#' 
#' @param l (numeric) A numeric number indicating the hyperparameter
#' (flexibility) of a specific kernel.
#' @param d (integer) For polynomial, d is the power; for matern, v = d + 1 /
#' 2; for rational, alpha = d.
#' @return \item{matrix_wise}{(function) A function calculating the relevance
#' of two matrices.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#' @export kernel_intercept
kernel_intercept <-
  function(l, d) {
    matrix_wise <- function(X1, X2, l, d) {
      1
    }
    matrix_wise
  }




#' Generating A Single Matrix-wise Function Using Linear
#' 
#' Generate matrix-wise functions for two matrices using linear kernel.
#' 
#' \bold{Polynomial Kernels} \deqn{k(x, x')=(x \cdot x')^d} We have intercept
#' kernel when \eqn{d=0}, and linear kernel when \eqn{d=1}.
#' 
#' @param l (numeric) A numeric number indicating the hyperparameter
#' (flexibility) of a specific kernel.
#' @param d (integer) For polynomial, d is the power; for matern, v = d + 1 /
#' 2; for rational, alpha = d.
#' @return \item{matrix_wise}{(function) A function calculating the relevance
#' of two matrices.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#' @export kernel_linear
kernel_linear <-
  function(l, d) {
    matrix_wise <- function(X1, X2, l, d) {
      X1 %*% t(X2)
    }
    matrix_wise
  }




#' Generating A Single Matrix-wise Function Using Polynomial
#' 
#' Generate matrix-wise functions for two matrices using polynomial kernel.
#' 
#' \bold{Polynomial Kernels} \deqn{k(x, x')=(x \cdot x')^d} We have intercept
#' kernel when \eqn{d=0}, and linear kernel when \eqn{d=1}.
#' 
#' @param l (numeric) A numeric number indicating the hyperparameter
#' (flexibility) of a specific kernel.
#' @param d (integer) For polynomial, d is the power; for matern, v = d + 1 /
#' 2; for rational, alpha = d.
#' @return \item{matrix_wise}{(function) A function calculating the relevance
#' of two matrices.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#' @export kernel_polynomial
kernel_polynomial <-
  function(l, d) {
    matrix_wise <- function(X1, X2, l, d) {
      (X1 %*% t(X2) + 1) ^ d
    }
    matrix_wise
  }




#' Generating A Single Matrix-wise Function Using RBF
#' 
#' Generate matrix-wise functions for two matrices using rbf kernel.
#' 
#' \bold{Gaussian RBF Kernels} \deqn{k_{SE}(r)=exp\Big(-\frac{r^2}{2l^2}\Big)}
#' 
#' @param l (numeric) A numeric number indicating the hyperparameter
#' (flexibility) of a specific kernel.
#' @param d (integer) For polynomial, d is the power; for matern, v = d + 1 /
#' 2; for rational, alpha = d.
#' @return \item{matrix_wise}{(function) A function calculating the relevance
#' of two matrices.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#' @export kernel_rbf
kernel_rbf <-
  function(l, d) {
    matrix_wise <- function(X1, X2, l, d) {
      exp(-square_dist(X1, X2, l = sqrt(2) * l))
    }
    matrix_wise
  }




#' Generating A Single Matrix-wise Function Using Matern
#' 
#' Generate matrix-wise functions for two matrices using matern kernel.
#' 
#' \bold{Matern Kernels}
#' \deqn{k_{Matern}(r)=\frac{2^{1-\nu}}{\Gamma(\nu)}\Big(\frac{\sqrt{2\nu
#' r}}{l}\Big)^\nu K_\nu \Big(\frac{\sqrt{2\nu r}}{l}\Big)}
#' 
#' @param l (numeric) A numeric number indicating the hyperparameter
#' (flexibility) of a specific kernel.
#' @param d (integer) For polynomial, d is the power; for matern, v = d + 1 /
#' 2; for rational, alpha = d.
#' @return \item{matrix_wise}{(function) A function calculating the relevance
#' of two matrices.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#' @export kernel_matern
kernel_matern <-
  function(l, d) {
    matrix_wise <- function(X1, X2, l, d){
      r <- sqrt(square_dist(X1, X2))
      v <- d + 1 / 2
      s <- 0
      for (i in 0:d) {
        s <- s + factorial(d + i) / (factorial(i) * factorial(d - i)) *
          (sqrt(8 * v) * r / l) ^ (d - i)
      }
      exp(-sqrt(2 * v) * r / l) * gamma(d + 1) / gamma(2 * d + 1) * s
    }
    matrix_wise
  }




#' Generating A Single Matrix-wise Function Using Rational Quadratic
#' 
#' Generate matrix-wise functions for two matrices using rational kernel.
#' 
#' \bold{Rational Quadratic Kernels} \deqn{k_{RQ}(r)=\Big(1+\frac{r^2}{2\alpha
#' l^2}\Big)^{-\alpha}}
#' 
#' @param l (numeric) A numeric number indicating the hyperparameter
#' (flexibility) of a specific kernel.
#' @param d (integer) For polynomial, d is the power; for matern, v = d + 1 /
#' 2; for rational, alpha = d.
#' @return \item{matrix_wise}{(function) A function calculating the relevance
#' of two matrices.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#' @export kernel_rational
kernel_rational <-
  function(l, d) {
    matrix_wise <- function(X1, X2, l, d){
      r <- sqrt(square_dist(X1, X2))
      (1 + r ^ 2 / (2 * d * l ^ 2)) ^ (-d)
    }
    matrix_wise
  }




#' Generating A Single Matrix-wise Function Using Neural Network
#' 
#' Generate matrix-wise functions for two matrices using neural network kernel.
#' 
#' \bold{Neural Network Kernels} \deqn{k_{NN}(x,
#' x')=\frac{2}{\pi}sin^{-1}\Big(\frac{2\tilde{x}^T
#' \tilde{x}'}{\sqrt{(1+2\tilde{x}^T \tilde{x})(1+2\tilde{x}'^T
#' \tilde{x}')}}\Big)}
#' 
#' @param l (numeric) A numeric number indicating the hyperparameter
#' (flexibility) of a specific kernel.
#' @param d (integer) For polynomial, d is the power; for matern, v = d + 1 /
#' 2; for rational, alpha = d.
#' @return \item{matrix_wise}{(function) A function calculating the relevance
#' of two matrices.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#' @export kernel_nn
kernel_nn <-
  function(l, d) {
    matrix_wise <- function(X1, X2, l, d){
      X1 <- cbind(1, X1)
      X2 <- cbind(1, X2)
      X1s <- apply(X1, 1, crossprod)
      X2s <- apply(X2, 1, crossprod)
      X1m <- matrix(X1s, nrow = length(X1s), ncol = nrow(X2), byrow = FALSE)
      X2m <- matrix(X2s, nrow = nrow(X1), ncol = length(X2s), byrow = TRUE)
      s <- 2 * (X1 %*% t(X2)) / (sqrt((1 + 2 * X1m) * (1 + 2 * X2m)))
      asin(s)
    }
    matrix_wise
  }




#' Computing Square Distance between Two Sets of Variables
#' 
#' Compute Squared Euclidean distance between two sets of variables with the
#' same dimension.
#' 
#' 
#' @param X1 (matrix, n1*p0) The first set of variables.
#' @param X2 (matrix, n2*p0) The second set of variables.
#' @param l (numeric) A numeric number indicating the hyperparameter
#' (flexibility) of a specific kernel.
#' @return \item{kern}{(function) A function calculating the relevance of two
#' matrices.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#' @export square_dist
square_dist <- function(X1, X2 = NULL, l = 1) {
  X1 <- X1 / l
  X1s <- apply(X1, 1, crossprod)
  if (is.null(X2)) {
    X2 <- X1
    X2s <- X1s
  } else {
    if (ncol(X1) != ncol(X2)) {
      stop("dimensions of X1 and X2 do not match!")
    }
    X2 <- X2 / l
    X2s <- apply(X2, 1, crossprod)
  }
  dist <- -2 * X1 %*% t(X2)
  X1m <- matrix(X1s, nrow = length(X1s), ncol = nrow(X2), byrow = FALSE)
  X2m <- matrix(X2s, nrow = nrow(X1), ncol = length(X2s), byrow = TRUE)
  dist + X1m + X2m
}

