% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/g.R
\name{g}
\alias{g}
\title{Gaussian mean estimation}
\usage{
g(x, s, rho = 0, K = 10, C = 5, tol = 1e-05, maxit = 100, verbose = FALSE)
}
\arguments{
\item{x}{Gaussian sequence}

\item{s}{standard deviation}

\item{rho}{regularization parameter, closer to 0 means less regularization}

\item{K}{number of grid points, in the interval [x_j - C s, x_j + C s], over which to search for the jth tuning parameter, where C is define below}

\item{C}{the value of the constant C in the interval above}

\item{tol}{error tolerance for convergence of the unbiased risk estimate}

\item{maxit}{maximum number of allowable iterations}

\item{verbose}{should the value of SURE be reported at each iteration}
}
\value{
\item{t_hat}{values of tuning parameters t}
\item{theta_hat}{estimated values of means of Gaussian sequence}
}
\description{
Estimates the mean vector of a primary homoscedastic sequence of independent Gaussian observations under squared error loss. The optimal decision rule is estimated by directly minimizing an unbiased estimate of its risk.
}
\examples{
\donttest{
## generate data
n = 250
set.seed(1)
theta = rnorm(n)
x = theta + rnorm(n)
## loss of MLE
mean((theta - x)^2)
## loss of estimator incorporating side information
mean((theta - g(x, 1)$theta_hat)^2)
}

}
