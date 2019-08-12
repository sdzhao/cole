#' Gaussian mean estimation
#'
#' Estimates the mean vector of a primary homoscedastic sequence of independent Gaussian observations under squared error loss. The optimal decision rule is estimated by directly minimizing an unbiased estimate of its risk.
#'
#'
#' @param x Gaussian sequence
#' @param s standard deviation
#' @param rho regularization parameter, closer to 0 means less regularization
#' @param K number of grid points, in the interval \[x_j - C s, x_j + C s\], over which to search for the jth tuning parameter, where C is define below
#' @param C the value of the constant C in the interval above
#' @param tol error tolerance for convergence of the unbiased risk estimate
#' @param maxit maximum number of allowable iterations
#' @param verbose should the value of SURE be reported at each iteration
#'
#' @return
#' \item{t_hat}{values of tuning parameters t}
#' \item{theta_hat}{estimated values of means of Gaussian sequence}
#'
#' @examples
#' \donttest{
#' ## generate data
#' n = 250
#' set.seed(1)
#' theta = rnorm(n)
#' x = theta + rnorm(n)
#' ## loss of MLE
#' mean((theta - x)^2)
#' ## loss of estimator incorporating side information
#' mean((theta - g(x, 1)$theta_hat)^2)
#' }
#' 
#' @useDynLib cole
#' @export

g = function(x, s, rho=0,
             K=10, C=5,
             tol=1e-5, maxit=100, verbose=FALSE) {
    ## error checking
    if (length(s) != 1) {
        stop("s must be scalar")
    }
    
    ## minimize sure
    t_hat = .Call("g_min_sure",
                   as.numeric(x), as.numeric(s), as.numeric(rho),
                   as.integer(K), as.numeric(C),
                   as.numeric(tol), as.integer(maxit), as.integer(2), as.integer(verbose))
    
    ## report results
    n = length(x)
    theta_hat = .Call("g_rule",
                      as.numeric(x), as.numeric(s),
                      as.numeric(t_hat), as.numeric(rho))

    return(list(t_hat = t_hat, theta_hat = theta_hat))
}


#' Separable rule for Gaussian mean estimation
#'
#' Given tuning parameter vector t, returns corresponding estimate for the mean vector of a homoscedastic sequence of independent Gaussian observations
#'
#'
#' @param x primary Gaussian sequence
#' @param s standard deviation of primary sequence
#' @param t tuning parameter vector t1
#' @param rho regularization parameter, closer to 0 means less regularization
#'
#' @return estimated values of means of primary Gaussian sequence
#'
#' @examples
#' ## generate data
#' n = 250
#' set.seed(1)
#' theta = rnorm(n)
#' x = theta + rnorm(n)
#' ## loss of MLE
#' mean((theta - x)^2)
#' ## loss of oracle separable estimator
#' mean((theta - g_rule(x, 1, theta))^2)
#' @useDynLib cole
#' @export

g_rule = function(x, s, t, rho=0) {
    ## error checking
    if(length(x) != length(t)){
        stop("x and t must be the same length")
    }
    
    ret = .Call("g_rule",
                as.numeric(x), as.numeric(s),
                as.numeric(t), as.numeric(rho))
    return(ret)
}
