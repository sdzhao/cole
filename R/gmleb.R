#' General maximum likelihood empirical Bayes estimator of Jiang and Zhang (2009)
#'
#' Estimates mean vector of a homoscedastic sequence of independent Gaussian observations using the nonparametric maximum likelihood procedure of Jiang and Zhang (2009)
#'
#'
#' @param x1 Gaussian sequence
#' @param s1 standard deviation of Gaussian sequence
#' @param init initial values for masses of discrete prior, default is discrete uniform over support points
#' @param grid grid of support points for discrete prior, default is M equally spaced points between min(x1) and max(x1)
#' @param M number of grid points
#' @param tol error tolerance of convergence of log likelihood
#' @param maxit maximum number of iterations
#'
#'
#' @return
#' estimated values of means of primary Gaussian sequence
#'
#' @examples
#' ## generate data
#' n = 250
#' set.seed(1)
#' theta1 = rnorm(n)
#' x1 = theta1 + rnorm(n)
#' ## loss of MLE
#' mean((theta1 - x1)^2)
#' ## loss of GMLEB
#' mean((theta1 - gmleb(x1, 1)$theta_hat)^2)
#' @useDynLib cole
#' @export

gmleb = function(x1, s1,
                  init = NULL, grid = NULL, M = 300,
                  tol = 1e-5, maxit = 1000) {

  ## error checking
  if (length(s1) != 1) {
    stop("s1 must be scalar")
  }

  if (is.null(init[1])) {
    init = rep(1 / M, M)
  }
  if (is.null(grid[1])) {
    grid = seq(min(x1), max(x1), length = M)
  }

  ## em algorithm
  A = exp(-0.5 * (outer(x1, grid, "-")^2 / s1^2)) / s1
  f = init
  err = Inf
  n = length(x1)

  for (it in 1:maxit) {
    oldf = f
    f = t(A) %*% (1 / (A %*% oldf)) * oldf / n
    oldll = sum(log(A %*% oldf))
    err = abs(sum(log(A %*% f)) - oldll) / abs(oldll)
    if (err <= tol) {
      break
    }
  }

  return(list(theta_hat = drop(A %*% (f * grid) / A %*% f)))
}
