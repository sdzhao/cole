#' Kernel Density f-Estimator
#'
#' An implementation of the f-estimator from "Nonparametric empirical Bayes and
#' compound decision approaches to estimation of a high-dimensional vector of
#' normal means" (Brown and Greenshtein, 2009) for estimating the mean of a
#' homoscedastic sequence of independent Gaussian observations under squared
#' error loss. The various small modifications and default parameters that are
#' suggested in their paper are implemented.
#'
#' @param x Gaussian sequence
#' @param s standard deviation
#' @param bw scalar bandwidth for Gaussian kernel density estimate
#' @param Cn scalar threshold for to clip estimates in low-density areas (Inf for no thresholding)
#'
#' @return
#' \item{theta_hat}{estimated values of means of Gaussian sequence}
#' @export
#' @importFrom stats dnorm
#'
#' @examples
#' theta = rnorm(250)
#' x = theta + rnorm(250)
#' res = f.kde(x, s = 1)
#' mean((theta - res$theta_hat)^2)
f.kde <- function(x, s, bw = s/sqrt(log(length(x))), Cn = s*sqrt(3*log(length(x)))) {
  # Gaussian kernel density estimates of f'(xi) / f(xi)
  # We can write as a matrix operation, but this formulation avoids declaring and inverting large matrices.
  R <- sapply(x, function(z) {
    f0 <- sum(stats::dnorm((x-z)/bw))
    f1 <- sum(stats::dnorm((x-z)/bw)*(x-z)/bw^2)
    f1/f0
  })

  # Threshold the ratio estimates and return
  list(theta_hat = x + (s^2+bw^2)*pmin(abs(R),Cn)*sign(R))
}
