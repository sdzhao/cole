#' Monotone Degree-Zero Weighted Trend Filter
#'
#' Estimate the mean of a homoscedastic sequence of independent Gaussian
#' observations under squared error loss. The optimal separable estimator is
#' estimated by minimizing a biased estimate of the risk under a monotone
#' non-decreasing constraint. Under the monotonicity constraint, the biased
#' risk estimate looks like a degree-zero weighted trend filter. The pooled
#' adjacent violators algorithm is used to fit the estimator. The default value
#' of `bw` is the optimal asymptotic rate (up to a constant).
#'
#' When `w` is "approx" the FFT transform is used to quickly generate weights
#' along `gr` in O(N log N). When `w` is "exact" a kernel density estimate is
#' used to generate weights along `gr` in O(N^2). Other weights can be specified
#' by passing a length `length(x)-1` numeric vector.
#'
#' If `knots = NULL` the knots will be placed at the mid-point of each
#' consecutive value of `x`: `sort(x)`. This is an approximation, but one that
#' works well in general and saves a lot of compute time. The correct knots are
#' placed at the minimum of the density between each consecutive observation.
#' Users can specify their own grid of knots by providing a length `length(x)-1`
#' numeric vector for `knots`.
#'
#' @param x Gaussian sequence
#' @param s standard deviation
#' @param bw scalar bandwidth for Gaussian kernel density estimate
#' @param w sequence of weights (see details for specification)
#' @param knots locations of knots
#' @param ... additional parameters passed to stats::density()
#'
#' @return
#' \item{theta_hat}{estimated values of means of Gaussian sequence}
#' \item{x}{original Gaussian sequence}
#' \item{s}{known standard deviation}
#' \item{bw}{bandwidth used for weights or NULL if user specified weights}
#' \item{w}{sequence of weights used}
#' \item{knots}{location of knots}
#' \item{risk.est}{estimated risk of the fitted estimator}
#' \item{...}{additional parameters passed to stats::density()}
#' @export
#' @importFrom stats density approx dnorm
#'
#' @examples
#' # basic usage
#' theta = rnorm(250)
#' x = theta + rnorm(250)
#' res = wtf0(x, s = 1)
#' mean((theta - res$theta_hat)^2)
#'
#' # fit visualization
#' plot(x, x, main = "WTF0") # observed
#' points(theta ~ x, col = "green") # true
#' points(theta_hat ~ x, as.data.frame(res[1:2]), col = "red") # estimated
#'
#' # alternate usage
#' ## w exact vs approx (very close)
#' all.equal(
#'   wtf0(x, s = 1, w = "exact"),
#'   wtf0(x, s = 1, w = "approx")
#' )
#'
#' ## Left knots
#' res2 = wtf0(x, s = 1, bw = 0.2, knots = sort(x)[-250]+min(diff(sort(x)))/100)
#' mean((theta - res2$theta_hat)^2)
#'
#' ## empirical weights
#' res3 = wtf0(x, s = 1, w = rep(1/250, 249))
#' mean((theta - res3$theta_hat)^2)
#'
#' ## extra parameters (better approximation and different bandwidth)
#' ## bw can only be a string when `w = "approx"` (see stats::density() for more
#' ## details).
#' res4 = wtf0(x, s = 1, w = "approx", n = 1048, bw = "SJ")
#' mean((theta - res4$theta_hat)^2)
wtf0 <- function(x, s, bw = s*length(x)^(-1/6), w = "approx", knots = NULL, ...) {
  # set-up
  n <- length(x)
  xs <- sort(x)

  # define grid of knots
  if (is.null(knots)) {
    knots <- xs[-n] + diff(xs)/2
  } else if (is.numeric(knots) & length(knots) == n-1) {
    knots <- sort(knots) # not needed, but helps de-bugging
  } else {
    stop("knots must be a numeric vector of values between sort(x).")
  }

  # define weights
  if (is.character(w)) {
    if (w == "approx") {
      w = stats::approx(stats::density(x = xs, bw = bw, kernel = "gaussian", ...), xout = knots)$y
    } else if (w == "exact") {
      w = colMeans(stats::dnorm(outer(xs, knots, FUN = "-"), sd = bw))
    } else {
      stop("w must be in c('approx', 'exact'), or be a vector of weights corresponding to knots between each value of sort(x).")
    }
  } else if (is.numeric(w) & length(knots) == n-1) {
    # values not used in estimation
    knots = NULL
    bw = NULL
  } else {
    stop("w must be in c('approx', 'exact'), or be a vector of weights corresponding to knots between each value of sort(x).")
  }

  # create the transformation/shift vector, absorbing the negative sign.
  # shift = t(D) %*% w (where D is the n-1 x n differences matrix in the trend filter penalty)
  shift <- n*(s^2)*c(w[1], diff(w), -w[n-1])

  # absorb linear penalty and use PAVA to fit monotone least-squares problem
  theta_hat <- pava_ls(xs + shift)

  # return
  list(
    theta_hat = theta_hat[rank(x, ties.method = "first")],
    x = x,
    s = s,
    bw = bw,
    w = w,
    knots = knots,
    risk.est = mean((xs - theta_hat)^2) + 2*(s^2)*sum(w * diff(theta_hat)) - (s^2),
    ...
  )
}




#' Optimal Knots for wtf0
#'
#' Compute the optimal knot placement for `wtf0()` by minimize f_h(t) between
#' each consecutive pair of observations. These knots are optimal for the risk
#' estimate, but may not produce a better minimizer of the true risk.
#'
#' @param x Gaussian sequence
#' @param bw scalar bandwidth for Gaussian kernel density estimate
#'
#' @return vector of optimal knots given x and bw (one element shorter than `x`)
#' @export
#' @importFrom stats optimize
#'
#' @examples
#' set.seed(1)
#' theta = rnorm(250)
#' x = theta + rnorm(250)
#' bw = 0.25
#'
#' # optimal knots
#' wtf0(x, 1, bw = bw, knots = optimal.knots(x, bw))$risk.est
#'
#' # approximate knots
#' wtf0(x, 1, bw = bw)$risk.est
optimal.knots <- function(x, bw) {
  x <- sort(x) # pre-sort xi

  sapply(seq_len(length(x) - 1), function(i) {
    stats::optimize(function(t) mean(dnorm(x - t, sd = bw)), interval = x[i:(i+1)])$minimum
  })
}
