#' Gaussian mean estimation with Gaussian side information
#'
#' Estimates the mean vector of a primary homoscedastic sequence of independent Gaussian observations under squared error loss by leveraging side information in the form of an auxiliary homoscedastic sequence of independent Gaussians.
#' 
#'
#' @param x1 primary Gaussian sequence
#' @param s1 standard deviation of primary sequence
#' @param x2 auxiliary Gaussian sequence of side information
#' @param s2 standard deviation of auxiliary sequence
#' @param rho regularization parameter, closer to 0 means less regularization
#' @param K number of grid points, in the interval [x_j1 - C s_1, x_j1 + C s1], over which to search for the jth tuning parameter, where C is define below
#' @param C the value of the constant C in the interval above
#' @param tol error tolerance for convergence of the unbiased risk estimate
#' @param maxit maximum number of allowable iterations
#'
#' @return
#' \item{t1_hat}{values of tuning parameters t1}
#' \item{t2_hat}{values of tuning parameters t2}
#' \item{theta_hat}{estimated values of means of primary Gaussian sequence}
#'
#' @examples
#' \donttest{
#' ## generate data
#' n = 250
#' set.seed(1)
#' theta1 = rnorm(n)
#' theta2 = theta1
#' x1 = theta1 + rnorm(n)
#' x2 = theta2 + rnorm(n)
#' ## loss of MLE
#' mean((theta1 - x1)^2)
#' ## loss of estimator incorporating side information
#' mean((theta1 - gg(x1, 1, x2, 1)$theta_hat)^2)
#' }
#'
#' @useDynLib cole
#' @export

gg = function(x1, s1, x2, s2, rho = 0,
              K = 10, C = 5,
              tol = 1e-5, maxit = 100) {
    ## error checking
    if (length(x1) != length(x2)) {
        stop("x1 and x2 must be the same length")
    }
    if (length(s1) != 1 || length(s2) != 1) {
        stop("s1 and s2 must be scalar")
    }
    if (rho < 0) {
        stop("rho must be positive")
    }
    
    ## minimize sure
    tmp = .Call("gg_min_sure",
                as.numeric(x1), as.numeric(s1),
                as.numeric(x2), as.numeric(s2),
                as.numeric(rho),
                as.integer(K), as.numeric(C),
                as.numeric(tol), as.integer(maxit))
    
    ## report results
    n = length(x1)
    t1_hat = tmp[1:n]
    t2_hat = tmp[(n + 1) : (2 * n)]
    theta_hat = .Call("gg_rule",
                      as.numeric(x1), as.numeric(s1),
                      as.numeric(x2), as.numeric(s2),
                      as.numeric(t1_hat), as.numeric(t2_hat),
                      as.numeric(rho))

    return(list(t1_hat = t1_hat, t2_hat = t2_hat, theta_hat = theta_hat))
}

#' Separable rule for Gaussian mean estimation with Gaussian side information
#'
#' Given tuning parameter vectors t1 and t2, returns corresponding estimate for the mean vector of a primary homoscedastic sequence of independent Gaussian observations that leverages side information in the form of an auxiliary homoscedastic sequence of independent Gaussians.
#' 
#'
#' @param x1 primary Gaussian sequence
#' @param s1 standard deviation of primary sequence
#' @param x2 auxiliary Gaussian sequence of side information
#' @param s2 standard deviation of auxiliary sequence
#' @param t1 tuning parameter vector t1
#' @param t2 tuning parameter vector t2
#' @param rho regularization parameter, closer to 0 means less regularization
#'
#' @return
#' \item{theta_hat}{estimated values of means of primary Gaussian sequence}
#'
#' @examples
#' ## generate data
#' n = 250
#' set.seed(1)
#' theta1 = rnorm(n)
#' theta2 = theta1
#' x1 = theta1 + rnorm(n)
#' x2 = theta2 + rnorm(n)
#' ## loss of MLE
#' mean((theta1 - x1)^2)
#' ## loss of oracle separable estimator
#' mean((theta1 - gg_rule(x1, 1, x2, 1, theta1, theta2)$theta_hat)^2)
#'
#' @useDynLib cole
#' @export

gg_rule = function(x1, s1, x2, s2, t1, t2, rho = 0) {
    ## error checking
    if (length(x1) != length(x2)) {
        stop("x1 and x2 must be the same length")
    }
    if (length(s1) != 1 || length(s2) != 1) {
        stop("s1 and s2 must be scalar")
    }
    if(length(x1) != length(t1) || length(x2) != length(t2)){
        stop("x1/x2 and t1/t2 must be the same length")
    }
    
    ret = .Call("gg_rule",
                as.numeric(x1), as.numeric(s1),
                as.numeric(x2), as.numeric(s2),
                as.numeric(t1), as.numeric(t2), as.numeric(rho))
    return(list(theta_hat = ret))
}
