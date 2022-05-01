#' Plug-in Estimate Oracle TV(d')
#'
#' A plug-in estimate of the total variation of the derivative of the oracle
#' decision rule, d*. This is calculated by plugging-in X_i for theta_i in
#' d*'', then R's numerical integration is used to estimate TV(d*').
#'
#' @param x Gaussian sequence
#' @param s standard deviation
#' @param ... other parameters passed into stats::integrate()
#'
#' @return an estimate of oracle TV(d')
#' @export
#' @importFrom stats dnorm integrate
#'
#' @examples
#' theta = rnorm(250)
#' TV1.oracle.plugin(theta, 1)
TV1.oracle.plugin <- function(x, s, ...) {
  abs.d2 <- function(y) {
    f0 <- rowMeans(stats::dnorm(outer(y, x, "-"), sd = s))
    f1 <- rowMeans(outer(y, x, function(a,b,sd) { (b-a)/sd/sd * stats::dnorm(a-b, sd = sd) }, sd = s))
    f2 <- rowMeans(outer(y, x, function(a,b,sd) { ((b-a)^2 - sd^2)/(sd^4) * stats::dnorm(a-b, sd = sd) }, sd = s))
    f3 <- rowMeans(outer(y, x, function(a,b,sd) { ((b-a)^3 - 3*sd^2*(b-a))/(sd^6) * stats::dnorm(a-b, sd = sd) }, sd = s))

    abs((s^2) * (f3 / f0) - 3 * (s^2) * (f2 * f1 / f0^2) + 2 * (s^2) * (f1 / f0)^3)
  }

  l <- min(x) - 3 * s * sqrt(2 * log(2 * length(x)))
  u <- max(x) + 3 * s * sqrt(2 * log(2 * length(x)))
  stats::integrate(abs.d2, lower = l, upper = u, ...)$value
}




#' High-probability Bound on Oracle TV(d')
#'
#' A high-probability upper bound for the total variation of the derivative of
#' the oracle decision rule, d*. This is calculated by first bounding TV(d*')
#' using the range of unobserved means, then estimating the range by expanding
#' the range of the observations.
#'
#' @param x Gaussian sequence
#' @param s standard deviation
#'
#' @return an upper bound on oracle TV(d')
#' @export
#'
#' @examples
#' theta = rnorm(250)
#' TV1.oracle.highprob(theta, 1)
TV1.oracle.highprob <- function(x, s) {
  s^(-2) * (max(x) - min(x) + 2*s*sqrt(3*log(length(x))))^2
}




#' Monotone Degree-One Weighted Trend Filter
#'
#' Estimate the mean of a homoscedastic sequence of independent Gaussian
#' observations under squared error loss. The optimal separable estimator is
#' estimated by minimizing an unbiased estimate of the risk under a monotone
#' non-decreasing constraint. Under the monotonicity constraint, the unbiased
#' risk estimate looks like a degree-one weighted trend filter.
#'
#' The use of the MOSEK backend or solver required a properly installed `Rmosek`
#' package. The latest installation instructions found here:
#' https://docs.mosek.com/latest/install/index.html and
#' https://docs.mosek.com/latest/rmosek/install-interface.html.
#'
#' The CVXR and MOSEK backends are designed for internal use only and do not
#' validate inputs.
#'
#' @param x Gaussian sequence
#' @param s standard deviation
#' @param tau a scalar bound on TV(d')
#' @param M number of evenly spaced knots (overrides knots if both specified)
#' @param knots locations of knots
#' @param backend one of "CVXR" or "MOSEK"
#' @param solver solver backend for CVXR
#' @param ... additional parameters passed to CVXR::psolve()
#' (if CVXR backend is used)
#'
#' @return
#' \item{theta_hat}{estimated values of means of Gaussian sequence}
#' \item{x}{original Gaussian sequence}
#' \item{s}{known standard deviation}
#' \item{SURE}{value of minimum risk estimate}
#' \item{tau}{the user specified bound on TV(d')}
#' \item{TV1}{the value of TV(d') for minimizer}
#' \item{intercept}{intercept of fitted estimator}
#' \item{slopes}{slopes of fitted estimator}
#' \item{knots}{location of knots}
#' \item{backend}{which version of wtf1 was used}
#' \item{solver}{solver used by CVXR}
#' \item{...}{additional inputed parameters}
#' @export
#'
#' @examples
#' # basic usage
#' set.seed(1)
#' theta = rnorm(250)
#' x = theta + rnorm(250)
#' res = wtf1(x, s = 1, tau = 1.2, knots = -2:2)
#' mean((theta - res$theta_hat)^2)
wtf1 <- function(x, s, tau = TV1.oracle.plugin(x = x, s = s), M = ceiling(sqrt(length(x))),
                 knots = NULL, backend = "CVXR", solver = "OSQP", ...) {
  # define knots
  ## sorting is needed for monotonicity constraint
  if (is.null(knots) & is.numeric(M)) {
    knots <- seq(from = min(x), to = max(x), length.out = M+2)[-c(1,M+2)]
  } else if (is.numeric(knots)) {
    knots <- sort(knots)
  } else {
    stop("knots not defined correctly, check values of knots and M.")
  }

  if (backend == "CVXR") {
    res <- wtf1.cvxr(x = x, s = s, tau = tau, knots = knots, solver = solver, ...)
  } else if (backend == "MOSEK") {
    res <- wtf1.mosek(x = x, s = s, tau = tau, knots = knots)
  } else {
    stop("backend should be one of 'CVXR' or 'MOSEK'.")
  }

  res
}




#' wtf1 CVXR backend
#'
#' @importFrom CVXR Variable Problem Minimize psolve
#'
#' @rdname wtf1
wtf1.cvxr <- function(x, s, tau, knots = NULL, solver = "OSQP", ...) {
  # define design matrix for the estimator and its derivative
  D <- outer(x, knots, function(a,b) pmax(a-b, 0))
  dD <- 1 * (D > 0)

  # define optimization problem keeping intercept separate from slopes
  beta0 <- CVXR::Variable()
  betas <- CVXR::Variable(length(knots))

  problem <- CVXR::Problem(
    objective = CVXR::Minimize(mean((x - beta0 - D %*% betas)^2) +
                                 2 * (s^2) * mean(dD %*% betas) - (s^2)),
    constraints = list(
      sum(abs(betas)) <= tau, # TV(d') <= tau
      cumsum(betas) >= 0, # monotone non-decreasing
      sum(betas) == 0 # saturation
    )
  )

  # solve constrained optimization problem
  res <- CVXR::psolve(problem, solver = solver, ...)

  # return estimated means and other problem and solver specifications
  list(
    theta_hat = res$getValue(beta0) + drop(D %*% res$getValue(betas)),
    SURE = res$value,
    tau = tau,
    TV1 = sum(abs(res$getValue(betas))),
    intercept = res$getValue(beta0),
    slopes = drop(res$getValue(betas)),
    knots = knots,
    backend = "CVXR",
    solver = res$solver,
    ...
  )
}



#' wtf1 MOSEK backend
#'
#' @rdname wtf1
wtf1.mosek <- function(x, s, tau, knots = NULL) {
  # Check for Rmosek installation
  if (!requireNamespace("Rmosek", quietly = TRUE)) {
    stop("Rmosek is required to use this function. Consider using the CVXR version instead.", call. = FALSE)
  }

  # set-up
  n <- length(x)
  p <- length(knots)
  D <- pmax(outer(x, knots, "-"), 0) # design matrix without intercept
  Ip <- diag(p)

  # define MOSEK optimization problem
  ## linear constraints in [a; b; u; t1]
  ## a: intercept; b: slopes; u: used for TV(d') bound; t1: used for conic
  prob <- list(sense = "min")
  ## linear constraint matrix, A
  ## A [a; b; u; t1] <,>,= b
  prob$A <- rbind(
    c(0, rep(1, p), rep(0, p), 0),
    cbind(0, lower.tri(Ip, diag = TRUE), matrix(0, p, p), 0),
    cbind(0, Ip, Ip, 0),
    cbind(0, -Ip, Ip, 0),
    c(0, rep(0, p), rep(1, p), 0)
  )
  ## constraints on A [a; b; u; t1]
  prob$bc <- rbind(
    c(0, rep(0, p), rep(0, p), rep(0, p), -Inf),
    c(0, rep(Inf, p), rep(Inf, p), rep(Inf, p), tau)
  )
  ## constraints on [a; b; u; t1]
  prob$bx <- rbind(
    rep(-Inf, 2*p+2),
    rep(Inf, 2*p+2)
  )
  ## linear min problem c'[a; b; u; t1]
  prob$c <- c(0, 2*(s^2)*colSums(D > 0), rep(0, p), 1)

  ## affine conic constraint
  ## F [a; b; u; t1] + g = [t1; 0.5; a + D b - x] in Q^{n+2}
  prob$F <- rbind(
    c(0, rep(0, p), rep(0, p), 1),
    c(0, rep(0, p), rep(0, p), 0),
    cbind(1, D, matrix(0, n, p), 0)
  )
  prob$g <- c(0, 0.5, -x)
  prob$cones <- matrix(list("RQUAD", n+2, NULL))
  rownames(prob$cones) <- c("type","dim","conepar")

  # solve constrained optimization problem
  res <- Rmosek::mosek(prob, list(verbose=1))

  # return estimated means and other problem and solver specifications
  list(
    theta_hat = res$sol$itr$xx[1] + drop(D %*% res$sol$itr$xx[2:(p+1)]),
    SURE = mean((x - res$sol$itr$xx[1] - D %*% res$sol$itr$xx[2:(p+1)])^2) +
      2*(s^2)*mean((D > 0) %*% res$sol$itr$xx[2:(p+1)]) - (s^2),
    tau = tau,
    TV1 = sum(abs(res$sol$itr$xx[2:(p+1)])),
    intercept = res$sol$itr$xx[1],
    slopes = res$sol$itr$xx[2:(p+1)],
    knots = knots,
    backend = "MOSEK",
    solver = "MOSEK"
  )
}




#' Monotone Degree-One Weighted Trend Filter Cross-Validation
#'
#' Use a fixed set of knots and a random partition of the data to cross-validate
#' the risk estimate for a provided list of tau values. If `nfolds = length(x)`
#' leave-one-out cross-validation is used. If `fit.best = TRUE` (default) then
#' a final model is fit using all of the data and best value of tau. A seed is
#' available to set for reproducible partitions.
#'
#' @param x Gaussian sequence
#' @param s standard deviation
#' @param tau.options a vector of candidate values for TV(d')
#' @param nfolds number of cross-validation folds
#' @param seed random seed used for reproducible folds
#' @param fit.best boolean (default TRUE) whether to fit all of `x` with the
#' best candidate value of TV(d')
#' @param ... additional parameters passed to wtf1() such as backend and knots
#'
#' @return
#' \item{theta_hat}{estimated values of means of Gaussian sequence}
#' \item{x}{original Gaussian sequence}
#' \item{s}{known standard deviation}
#' \item{SURE}{value of minimum risk estimate}
#' \item{tau}{the user specified bound on TV(d')}
#' \item{TV1}{the value of TV(d') for minimizer}
#' \item{intercept}{intercept of fitted estimator}
#' \item{slopes}{slopes of fitted estimator}
#' \item{knots}{location of knots}
#' \item{backend}{which version of wtf1 was used}
#' \item{solver}{solver used by CVXR}
#' \item{...}{additional inputed parameters}
#' \item{cv.risks}{data.frame containing cross-validation parameters and risks}
#' \item{nfold}{number of folds for cross-validation}
#' \item{seed}{seed for reproducible fold partitions}
#' @export
#'
#' @examples
#' # basic usage
#' set.seed(1)
#' theta = rnorm(250)
#' x = theta + rnorm(250)
#' res = cv.wtf1(x, s = 1, tau.options = 1:2, knots = -2:2)
#' mean((theta - res$theta_hat)^2)
cv.wtf1 <- function(x, s, tau.options, nfolds = 5, seed = NULL, fit.best = TRUE, ...) {
  # check nfolds
  if (nfolds < 1 | nfolds > length(x)) {
    stop("nfolds should be bewteen 1 and length(x)")
  }
  if (length(x) %% nfolds != 0) {
    warning("length(x) cannot be evenly devided into nfolds partitions")
  }

  # run nfold cross-validation for each tau with a random folds
  set.seed(seed = seed)
  folds <- split(x = sample(seq_along(x)), f = seq_along(x) %% nfolds)
  cv.risk <- sapply(tau.options, function(tau) {
    mean(sapply(folds, function(i) {
      fit <- wtf1(x = x[-i], s = s, tau = tau, ...)
      pred.x <- fit$intercept + pmax(outer(x[i], fit$knots, "-"), 0) %*% fit$slopes
      (x[i] - pred.x)^2 + 2*(s^2)*drop(crossprod(fit$slopes, outer(fit$knots, x[i], FUN = "<"))) - (s^2)
    }))
  })

  # fit best tau if requested
  if (fit.best) {
    res <- wtf1(x = x, s = s, tau = tau.options[which.min(cv.risk)], ...)
  } else {
    res <- list()
  }

  # return results
  res$cv.risks <- data.frame(tau = tau.options, cv.risk = cv.risk)
  res$nfolds <- nfolds
  res$seed <- seed
  res
}
