# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Compound decision method for two-sieded models
#'
#' EM algorithm for predicting unseen data for two-sided linear model based on candidate parameters estimated by simple models
#' 
#'
#' @param y outcomes for training with dimension n by q
#' @param x features for training with dimension n by p
#' @param ss condidates sigma square estimated by other mithods (lasso, ridge, camel)
#' @param B condidates beta estimated by other mithods (lasso, ridge, camel)
#' @param tol error tolerance for convergence
#' @param maxit maximum number of allowable iterations
#' @param p number of features
#' @param q number of outcomes
#' @param n number of observations for training
#' @param newx new observation for prediction
#'
#' @return
#' \item{f}{estimation of prior of beta}
#' \item{A}{likelihood of training data}
#' \item{bs}{candidates beta}
#' \item{esty}{estimation for based on newx}
#'
#' @examples
#' \donttest{
#' p = 10
#' q = 5
#' n = 50
#' 
#' x = matrix(rnorm(n*p,0,10), n, p)
#' beta = matrix(rnorm(p*q,0,10), q, p)
#' e = matrix(rnorm(n*q,0,0.1),n,q)
#' y = x %*% t(beta) + e
#' s2 = matrix(rep(0.1,q), q, 1)
#' tol = 0.001
#' maxit = 1000
#' x_test = matrix(rnorm(n*p,0,1), n, p)
#' 
#' output = comte(y, x, s2, beta, tol, maxit, p, q, n, x_test)
#' 
#' }
#'
#' @useDynLib cole
#' @export

comte <- function(y, x, ss, B, tol, maxit, p, q, n, newx) {
    .Call('_cole_comte', PACKAGE = 'cole', y, x, ss, B, tol, maxit, p, q, n, newx)
}


#' Compound decision method for two-sieded models with minimal sigma square
#'
#' EM algorithm for predicting unseen data for two-sided linear model based on candidate parameters estimated by simple models. Sigma square is taken to be the minimum for better accuracy in some situations
#' 
#'
#' @param y outcomes for training with dimension n by q
#' @param x features for training with dimension n by p
#' @param ss condidates sigma square estimated by other mithods (lasso, ridge, camel)
#' @param B condidates beta estimated by other mithods (lasso, ridge, camel)
#' @param tol error tolerance for convergence
#' @param maxit maximum number of allowable iterations
#' @param p number of features
#' @param q number of outcomes
#' @param n number of observations for training
#' @param newx new observation for prediction
#'
#' @return
#' \item{f}{estimation of prior of beta}
#' \item{A}{likelihood of training data}
#' \item{bs}{candidates beta}
#' \item{esty}{estimation for based on newx}
#'
#' @examples
#' \donttest{
#' p = 10
#' q = 5
#' n = 50
#' 
#' x = matrix(rnorm(n*p,0,10), n, p)
#' beta = matrix(rnorm(p*q,0,10), q, p)
#' e = matrix(rnorm(n*q,0,0.1),n,q)
#' y = x %*% t(beta) + e
#' s2 = matrix(rep(0.1,q), q, 1)
#' tol = 0.001
#' maxit = 1000
#' x_test = matrix(rnorm(n*p,0,1), n, p)
#' 
#' output = comte_min(y, x, s2, beta, tol, maxit, p, q, n, x_test)
#' 
#' }
#'
#' @useDynLib cole
#' @export

comte_min <- function(y, x, ss, B, tol, maxit, p, q, n, newx) {
    .Call('_cole_comte_min', PACKAGE = 'cole', y, x, ss, B, tol, maxit, p, q, n, newx)
}
