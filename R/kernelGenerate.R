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
#' @param method A character string indicating which kernel is to be computed.
#' @param Sigma The covariance matrix for neural network kernel.
#' @param l A numeric number indicating the hyperparameter (flexibility) of a
#' specific kernel.
#' @param p For polynomial, p is the power; for matern, v = p + 1 / 2; for
#' rational, alpha = p.
#' @return \item{Kern}{A function indicating the generated kernel.}
#' @author Wenying Deng
#' @references The MIT Press. Gaussian Processes for Machine Learning, 2006.
#' @examples
#'
#'
#' ##kernelGenerate(method = "rbf", Sigma = 0, l = 1, p = 2)
#'
#' ##Kernlist <- NULL
#' ##Kernlist <- c(Kernlist, kernelGenerate('rbf', l = .6))
#' ##Kernlist <- c(Kernlist, kernelGenerate('rbf', l = 1))
#' ##Kernlist <- c(Kernlist, kernelGenerate('rbf', l = 2))
#'
#'
#' @export kernelGenerate
#'
#' @import mvtnorm MASS psych limSolve stats
kernelGenerate <-
  function(method = "rbf", Sigma = 0, l = 1, p = 2){

    if (method == "intercept")
      SE <- function(xp, xq, Sigma, l, p) 1
    else if (method == "linear")
      SE <- function(xp, xq, Sigma, l, p) t(xp) %*% xq
    else if (method == "polynomial")
      SE <- function(xp, xq, Sigma, l, p) (t(xp) %*% xq + 1) ^ p
    else if (method == "rbf")
      SE <- function(xp, xq, Sigma, l, p) exp(- sum((xp - xq) ^ 2) / (2 * l ^ 2))
    else if (method == "matern")
      SE <- function(xp, xq, Sigma, l, p){
        r <- sqrt(sum((xp - xq) ^ 2))
        v <- p + 1 / 2
        s <- 0
        for (i in 0:p) {
          s <- s + factorial(p + i) / (factorial(i) * factorial(p - i)) *
            (sqrt(8 * v) * r / l) ^ (p - i)
        }
        k_v <- exp(-sqrt(2 * v) * r / l) * gamma(p + 1) / gamma(2 * p + 1) * s
        return(k_v)
      }
    else if (method == "rational")
      SE <- function(xp, xq, Sigma, l, p){
        r <- sqrt(sum((xp - xq) ^ 2))
        k_a <- (1 + r ^ 2 / (2 * p * l ^ 2)) ^ (- p)
        return(k_a)
      }
    else if (method == "nn")
      SE <- function(xp, xq, Sigma, l, p){
        xp <- c(1, xp)
        xq <- c(1, xq)
        if(Sigma == 0) Sigma <- diag(length(xp))
        s <- 2 * t(xp) %*% Sigma %*% xq / (sqrt((1 + 2 * t(xp) %*% Sigma %*% xp)
                                                * (1 + 2 * t(xq) %*% Sigma %*% xq)))
        k_n <- asin(s)
        return(k_n)
      }
    else
      stop("method must be one kind of kernel!")

    Kern <- function(X2, X1) apply(X1, 1, function(xp){
      apply(X2, 1, function(xq){
        SE(xp, xq, Sigma, l, p)
      })
    })

    return(Kern)
  }
