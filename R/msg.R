#' A Compound Decision Approach to Covariance Matrix Estimation
#'
#' Estimates covariance matrix for data from Gaussian distribution
#'
#'
#' @param X n x p data matrix, each column is a feature
#' @param K number of clusters
#' @param centered whether to centralize each column of X
#' @param maxit maximum number of iterations
#' @param tol error tolerance of convergence of log likelihood
#' @param verbose whether to output the error in each iteration
#'
#'
#' @return
#' estimated p x p covariance matrix
#'
#' @examples
#' ## generate data
#' n = 100; p = 30
#' set.seed(1)
#' sigma.x = matrix(0.9,p,p)
#' diag(sigma.x) = 1
#' X = MASS::mvrnorm(n,rep(0,p),sigma.x)
#' ## loss of msg
#' norm(sigma.x - msg(X,K=p),'f')
#' @useDynLib cole
#' @export
#' 
msg = function(X, K=10, n_start=20, centered=FALSE, maxit = 200, tol = 1e-4, verbose = FALSE) {
  lp = ncol(X); ln=nrow(X)
  if (centered) {
    ln = ln-1
    X = scale(X, scale = FALSE)
  }
  S = t(X) %*% X 
  ## make grid
  sds = sqrt(diag(S / ln))
  row_sdmat = replicate(lp,sds)
  col_sdmat = t(row_sdmat)
  cor_mat = S / ln / row_sdmat / col_sdmat
  grid = matrix(nrow=lp*(lp-1)/2,ncol=3)
  grid[,1] = row_sdmat[lower.tri(row_sdmat, diag=F)]
  grid[,2] = col_sdmat[lower.tri(col_sdmat, diag=F)]
  grid[,3] = cor_mat[lower.tri(cor_mat, diag=F)]
  nd_grid = grid
  n_dis = nrow(unique(nd_grid))
  K1 = min(n_dis,K)
  
  # grid clusters for non-diagonal entries
  if (n_dis < K)
  {
    grid.l = as.matrix(unique(nd_grid))
  }
  else
  {
    cnd_grid = scale(nd_grid)
    ct_nd_grid = attr(cnd_grid, 'scaled:center')
    sc_nd_grid = attr(cnd_grid, 'scaled:scale')
    # set.seed(233)
    cl_nd = kmeans(cnd_grid, centers=K, nstart = n_start,iter.max = 100)
    cnd_centers = cl_nd$centers
    nd_centers = apply(cnd_centers, 1,function(x) {x * sc_nd_grid + ct_nd_grid})
    grid.l = t(nd_centers)
  }
  
  grid.u = matrix(NA, nrow=K1, ncol=3)
  grid.u[,1] = grid.l[,2]
  grid.u[,2] = grid.l[,1]
  grid.u[,3] = grid.l[,3]
  grid = rbind(grid.l, grid.u)
  
  cmbn = combn(lp, 2)
  cmbn = cmbn[c(2,1),]
  W.all = matrix(NA,nrow=ncol(cmbn),ncol = 4)
  W.all[,1] = S[cbind(cmbn[1,],cmbn[1,])]
  W.all[,2] = S[cbind(cmbn[1,],cmbn[2,])]
  W.all[,3] = W.all[,2]
  W.all[,4] = S[cbind(cmbn[2,],cmbn[2,])]
  det.W = W.all[,1] * W.all[,4] - W.all[,2]^2
  
  log.det.w = log(det.W)
  log.I1 = (ln-3)/2*log.det.w
  log.I4 = ln*log(2)+log(pi)/2+log(gamma(ln/2)) + log(gamma((ln-1)/2))
  
  D = apply(grid, 1, function(x) {
    ## off-diagonal
    sigma = matrix(c(x[1]^2, prod(x), prod(x), x[2]^2),nrow = 2)
    inv.sigma = solve(sigma)
    
    log.I2 = -(W.all %*% as.vector(inv.sigma)) / 2
    log.I3 = ln / 2 * log(det(sigma))
    tmp = log.I1 + log.I2 - log.I3 - log.I4
    off = exp(tmp)
    ## on-diagonal is scaled chi-square
    on = dchisq(diag(S)/x[1]^2, ln) * ln / x[1]^2
    return(c(off, on))
  })
  g = npmle_sym_mat(D, maxit, tol, verbose)
  ind = lp * (lp-1) / 2 
  
  ## estimated covariance matrix
  denom = D %*% g
  tmp_nd = D[1:ind,] %*% (g * apply(grid, 1, prod)) / denom[1:ind]
  tmp_d = D[(ind+1):(ind+lp),] %*% (g * grid[,1]^2) / denom[(ind+1):(ind+lp)]
  est_cov = matrix(0, lp, lp)
  est_cov[cbind(cmbn[1,],cmbn[2,])] = as.vector(tmp_nd)
  est_cov[cbind(cmbn[2,],cmbn[1,])] = as.vector(tmp_nd)
  diag(est_cov) = as.vector(tmp_d)
  
  return(est_cov)
}

#' npmle_sym_mat
#'
#' estimate the weights of grid points maximizing likelihood with given density matrix 
#'
#' @param D density matrix
#' @param maxit maximum number of iterations
#' @param tol error tolerance of convergence of log likelihood
#' @param verbose whether to output the error in each iteration
#'
#'
#' @return
#' estimated weights of grids
#' 
#' @useDynLib cole
#' @export
#' 
#' 
npmle_sym_mat = function (D, maxit = 200, tol = 1e-04, verbose = FALSE)  {
  lp = nrow(D)
  d = ncol(D)
  g = rep(1/d, d)
  D.l = D[,1:(d/2)]
  D.u = D[,(d/2+1):d]
  D.half = (D.l + D.u) / 2
  D.avg = cbind(D.half, D.half)
  oldll = 0
  for (it in 1:maxit) {
    tmp = D %*% g
    newll = sum(log(tmp))
    err = abs(newll - oldll)/abs(oldll)
    if (err <= tol) {
      break
    }
    else {
      oldll = newll
      g = crossprod(D.avg, 1/tmp) * g/lp
      if (verbose)
        cat("Iter", it, ", error", err, "\n")
    }
  }
  if (it == maxit) {
    warning("Maximum number of iterations reached.")
  }
  return(g)
}

