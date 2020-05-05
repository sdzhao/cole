#' Positive definiteness correction for a given symmetric sqaure matrix
#' 
#' This function is used to find an optimal positive definite matrix for a given symmetric matrix. First,
#' it proposes a d-length sequence of exponent parameters, which constructs the candidates of smallest 
#' positive eigenvalues with the minimal positive eigenvalue of the given matrix. The smallest positive 
#' eigenvalues is determined by a penalized minimal Frobenius norm equation. This method is implemented 
#' from Huang & Chao's paper (2017).
#' 
#' Reference: Huang, Chao, Daniel Farewell, and Jianxin Pan. "A calibration method for non-positive 
#' definite covariance matrix in multivariate data analysis." Journal of Multivariate Analysis 157 
#' (2017): 45-52. 
#' 
#' @param Sigma    matrix for calibration
#' @param d        number of proposed smallest positive eigenvalues
#' 
#' @return         calibrated positive definite matrix
#' 
#' @examples 
#' p=10
#' set.seed(1)
#' mat = matrix(runif(p^2),p,p)
#' mat = mat + t(mat)
#' diag(mat) = 4
#' print(eigen(mat)$values)
#' mat1 = posdef.correction(mat,d=20)
#' print(eigen(mat1)$values)
#' 
#' @export

posdef.correction = function(Sigma, d)
{
  p = nrow(Sigma)
  alphas = seq(0,10,length.out = d) # d-length sequence of exponent parameter 
  eig = eigen(Sigma) 
  sigma.eig = eig$values
  U = eig$vectors
  lam.min = min(sigma.eig[sigma.eig>0]) # find the minimal positive eigenvalue
  c = 10^(-alphas) * lam.min 
  loss = rep(0,d)
  for(i in 1:d)
  {
    lam = pmax(sigma.eig, c[i])
    Sigma1 = U %*% diag(lam) %*% t(U)
    loss[i] = norm(Sigma1-Sigma,'f')/sqrt(p)+alphas[i]
    
  }
  i.m = which.min(loss)
  lam = pmax(sigma.eig, c[i.m])
  Sigma1 = U %*% diag(lam) %*% t(U)
  return(Sigma1)
}
