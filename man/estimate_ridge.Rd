% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation.R
\name{estimate_ridge}
\alias{estimate_ridge}
\title{Estimating a Single Model}
\usage{
estimate_ridge(Y, X, K, lambda, compute_kernel_terms = TRUE,
  converge_thres = 1e-04)
}
\arguments{
\item{Y}{(vector of length n) Reponses of the dataframe.}

\item{X}{(dataframe, n*p) Fixed effects variables in the dataframe (could
contains several subfactors).}

\item{K}{(list of matrices) A nested list of kernel term matrices, 
corresponding to each kernel term specified in the formula for 
a base kernel function in kern_func_list.}

\item{lambda}{(numeric) A numeric string specifying the range of tuning parameter 
to be chosen. The lower limit of lambda must be above 0.}

\item{compute_kernel_terms}{(logic) Whether to computing effect for each individual terms.
If FALSE then only compute the overall effect.}

\item{converge_thres}{(numeric) The convergence threshold for computing kernel terms.}
}
\value{
\item{lambda}{(numeric) The selected tuning parameter based on the
estimated ensemble kernel matrix.}

\item{beta}{(matrix, p*1) Fixed effects estimator of the model.}

\item{alpha}{(matrix, n*length(K)) Random effects estimator for each kernel 
term specified in the formula.}

\item{proj_matrix}{(list of length 4) Estimated projection matrices of the
model.}
}
\description{
Estimating projection matrices and parameter estimates for a single model.
}
\details{
For a single model, we can calculate the output of gaussian process
regression, the solution is given by \deqn{\hat{\beta}=[X^T(K+\lambda
I)^{-1}X]^{-1}X^T(K+\lambda I)^{-1}y} \deqn{\hat{\alpha}=(K+\lambda
I)^{-1}(y-\hat{\beta}X)}.
}
\examples{



estimate_ridge(Y = CVEK:::model_matrices$y, 
X = CVEK:::model_matrices$X, K = CVEK:::K_ens, 
lambda = CVEK:::lambda_ens)



}
\references{
Andreas Buja, Trevor Hastie, and Robert Tibshirani. (1989) 
Linear Smoothers and Additive Models. Ann. Statist. Volume 17, Number 2, 453-510.
}
\author{
Wenying Deng
}
