#' Unconstrained Coordinate Linear Shrinkage method (COOLISH) 
#'
#' Linearly shrinks coordinate of least squared estimates to solve multivariate linear prediction problems.
#'
#' @param x_tr training predictor matrix (n_tr x p)
#' @param y_tr training response matrix (n_tr x q)
#' @param x_te testing predictor matrix (n_te x p)
#'
#' @return
#' \item{y}{estimated response matrix of testing set (n_te x q)}
#'
#' @examples
#' n_tr = 100
#' n_te = 10
#' p = 30
#' q = 1000
#' 
#' set.seed(1)
#' B = matrix(rnorm(p*q), p, q)
#' B[6:p, ] = 0
#' x_tr = matrix(rnorm(n_tr*p), n_tr, p)
#' x_te = matrix(rnorm(n_te*p), n_te, p)
#' 
#' y_tr = x_tr %*% B + matrix(rnorm(n_tr*q), n_tr, q)
#' y_te = x_te %*% B + matrix(rnorm(n_te*q), n_te, q)
#' 
#' y_coolish = coolish_uncon(x_tr, y_tr, x_te)
#' fit_ols = lm(y_tr~x_tr)
#' y_ols = cbind(1, x_te) %*% fit_ols$coefficients
#' 
#' mean((y_te-y_coolish)^2)
#' mean((y_te-y_ols)^2)
#' @export

coolish_uncon = function(x_tr,
                         y_tr,
                         x_te){
  
  n_te = nrow(x_te)
  q = ncol(y_tr)
  p = ncol(x_tr)
  y_hat = matrix(NA, n_te, q)
  
  x = cbind(1, x_tr)
  fit = fit = lm(y_tr ~ x_tr)
  V = solve(crossprod(x, x)) * sum(sapply(summary(fit), function(x) x$sigma^2))
  x_new = cbind(1, x_te)
  y_ols = x_new %*% coef(fit)
  for (j in 1:nrow(x_new)) {
    index_zero = as.numeric(which(x_new[j, ] == 0) + 1) 
    index_nonzero = setdiff(1:(p+2), index_zero)
    Z = cbind(1, t(coef(fit) * x_new[j, ]))
    L = c(0, (x_new[j,] %*% V) * x_new[j,])
    
    if(length(index_zero) >= 1){
      Z = Z[, index_nonzero]
      L = L[index_nonzero]
    }
    a = solve(crossprod(Z, Z)) %*% (crossprod(Z, y_ols[j, ]) - L)
    y_hat[j,] = Z %*% a 
  }
  return(y_hat)
}

#' Constrained Coordinate Linear Shrinkage method (COOLISH) 
#'
#' Linearly shrinks coordinate of least squared estimates with constrains to solve multivariate linear prediction problems.
#'
#' @param x_tr training predictor matrix (n_tr x p)
#' @param y_tr training response matrix (n_tr x q)
#' @param x_te testing predictor matrix (n_te x p)
#'
#' @return
#' \item{y}{estimated response matrix of testing set (n_te x q)}
#'
#' @examples
#' n_tr = 100
#' n_te = 10
#' p = 30
#' q = 1000
#' 
#' set.seed(1)
#' B = matrix(rnorm(p*q), p, q)
#' B[6:p, ] = 0
#' x_tr = matrix(rnorm(n_tr*p), n_tr, p)
#' x_te = matrix(rnorm(n_te*p), n_te, p)
#' 
#' y_tr = x_tr %*% B + matrix(rnorm(n_tr*q), n_tr, q)
#' y_te = x_te %*% B + matrix(rnorm(n_te*q), n_te, q)
#' 
#' y_coolish = coolish_con(x_tr, y_tr, x_te)
#' fit_ols = lm(y_tr~x_tr)
#' y_ols = cbind(1, x_te) %*% fit_ols$coefficients
#' 
#' mean((y_te-y_coolish)^2)
#' mean((y_te-y_ols)^2)
#' @import quadprog
#' @export

coolish_con = function(x_tr,
                       y_tr,
                       x_te){
  
  n_te = nrow(x_te)
  q = ncol(y_tr)
  p = ncol(x_tr)
  y_hat = matrix(NA, n_te, q)
  
  x = cbind(1, x_tr)
  fit = fit = lm(y_tr ~ x_tr)
  V = solve(crossprod(x, x)) * sum(sapply(summary(fit), function(x) x$sigma^2))
  x_new = cbind(1, x_te)
  y_ols = x_new %*% coef(fit)
  
  for (j in 1:nrow(x_new)) {
    index_zero = as.numeric(which(x_new[j, ] == 0) + 1)
    index_nonzero = setdiff(1:(p+2), index_zero)
    Z = cbind(1, t(coef(fit) * x_new[j, ]))
    L = c(0, (x_new[j,] %*% V) * x_new[j,])
    
    if(length(index_zero) >= 1){
      Z = Z[, index_nonzero]
      L = L[index_nonzero]
    }
    Z_p = ncol(Z) - 2
    Dmat = 2 * crossprod(Z) /q
    dvec = as.numeric(2 * (y_ols[j,] %*% Z - L)/q)
    Amat = diag(1, Z_p+2)[, 2:(Z_p+2)]
    bvec = rep(0, Z_p+1)
    
    qpres = solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec)
    a = qpres$solution
    y_hat[j,] = Z %*% a
  }
  
  return(y_hat)
}
