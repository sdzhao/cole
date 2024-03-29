% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R, R/msg.R
\name{msg}
\alias{msg}
\title{Matrix shrinkage by g-modeling method}
\usage{
msg(
  X,
  K = 10,
  n_start = 20,
  centered = FALSE,
  maxit = 200,
  tol = 1e-04,
  verbose = FALSE
)

msg(
  X,
  K = 10,
  n_start = 20,
  centered = FALSE,
  maxit = 200,
  tol = 1e-04,
  verbose = FALSE
)
}
\arguments{
\item{X}{n x p data matrix, each column is a feature}

\item{K}{number of clusters}

\item{centered}{whether to centralize each column of X}

\item{maxit}{maximum number of iterations}

\item{tol}{error tolerance of convergence of log likelihood}

\item{verbose}{whether to output the error in each iteration}

\item{d}{interger. It controls the number of support points. Default value is 10.}
}
\value{
\itemize{
\item \code{est_cov}: \eqn{p x p} matrix of estiamted covariance matrix.
}

estimated p x p covariance matrix
}
\description{
The function use shrinkage method to estimate covariance matrix with given multinormal data.
This method applies empirical method and approximate the optimal separable decision rule by
EM algorithm.

Estimates covariance matrix for data from Gaussian distributino
}
\examples{
\donttest{
## set parameters
p = 10
d = 20
n = 100
## generate multivariate normal data X, with mean of zeros and identity matrix as covariance matrix.
X = matrix(rnorm(n*p,0,1), n, p)
## estimate covariance with given data X
cov = msg(X=X, d=d)
}

## generate data
n = 100; p = 30
set.seed(1)
sigma.x = matrix(0.9,p,p)
diag(sigma.x) = 1
X = mvrnorm(n,sigma.x)
## loss of msg
norm(sigma.x - msg(X,K=p),'f')
}
