% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tuning.R
\name{tuning_GCVc}
\alias{tuning_GCVc}
\title{Calculating Tuning Parameters Using GCVc}
\usage{
tuning_GCVc(Y, X, K_mat, lambda)
}
\arguments{
\item{Y}{(vector of length n) Reponses of the dataframe.}

\item{X}{(dataframe, n*p) Fixed effects variables in the dataframe (could
contains several subfactors).}

\item{K_mat}{(matrix, n*n) Estimated ensemble kernel matrix.}

\item{lambda}{(numeric) A numeric string specifying the range of tuning parameter 
to be chosen. The lower limit of lambda must be above 0.}
}
\value{
\item{lambda0}{(numeric) The estimated tuning parameter.}
}
\description{
Calculate tuning parameters based on GCVc.
}
\details{
\bold{Generalized Cross Validation (small sample size)}

\deqn{\lambda_{GCVc}=\underset{\lambda \in \Lambda}{argmin}\Big\{log\;
y^{\star
T}(I-A_\lambda)^2y^\star-2log[1-\frac{tr(A_\lambda)}{n}-\frac{2}{n}]_+\Big\}}
}
\references{
Philip S. Boonstra, Bhramar Mukherjee, and Jeremy M. G. Taylor.
A Small-Sample Choice of the Tuning Parameter in Ridge Regression. July
2015.

Trevor Hastie, Robert Tibshirani, and Jerome Friedman. The Elements of
Statistical Learning: Data Mining, Inference, and Prediction, Second
Edition. Springer Series in Statistics. Springer- Verlag, New York, 2
edition, 2009.

Hirotogu Akaike. Information Theory and an Extension of the Maximum
Likelihood Princi- ple. In Selected Papers of Hirotugu Akaike, Springer
Series in Statistics, pages 199–213. Springer, New York, NY, 1998.

Clifford M. Hurvich and Chih-Ling Tsai. Regression and time series model
selection in small samples. June 1989.

Hurvich Clifford M., Simonoff Jeffrey S., and Tsai Chih-Ling. Smoothing
parameter selection in nonparametric regression using an improved Akaike
information criterion. January 2002.
}
\author{
Wenying Deng
}
