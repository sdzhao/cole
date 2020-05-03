#include <Rmath.h>
#include <cmath>
#include <stdio.h>      
#include <RcppArmadillo.h>
#include <time.h>

using namespace Rcpp;
using namespace arma;

// same as rep(x,times) function in R, repeat a vector x for `times` times 
// [[Rcpp::export]]
arma::vec rep_time(arma::vec s, int times)
{
  int l = s.size();
  int total_l = l * times;
  arma::vec result(total_l);
  for(int i=0;i<times;i++)
  {
    result.subvec(i*l,(i+1)*l-1) = s;
  }
  return result;
}

// same as rep(x,each) function in R, repeat each element in vector x for `each` times 
// [[Rcpp::export]]
arma::vec rep_each(arma::vec s, int each)
{
  int l = s.size();
  int total_l = l * each;
  arma::vec tmp(each);
  arma::vec result(total_l);
  for(int i=0;i<l;i++)
  {
    result.subvec(i*each,(i+1)*each-1) = tmp.fill(s[i]);
  }
  return result;
}

// EM algorithm for nonparametric maximam likelihood estimation 
// [[Rcpp::export]]
arma::vec npmle_c(arma::mat D,int maxit=200, double tol=0.0001,bool verbose=false)
{
  int p = D.n_rows;
  int d = D.n_cols;
  arma::vec g(d);
  g.fill((double)(1.0/d));
  double oldll=0;
  double newll;
  double err;
  arma::vec tmp(p);
  for(int it=0;it<maxit;it++)
  {
    tmp = D * g;
    newll = accu(log(tmp));
    if(abs(oldll)>0){
      err= std::abs(newll - oldll) / std::abs(oldll);
    }
    else
    {
      err = R_PosInf;
    }
    if(err <= tol)
    {
      break;
    }
    else{
      oldll = newll;
      g = ((1/tmp).t() * D).t() % g / p;
      if(verbose)
      {
        printf("Iter %d, error %f \n", it, err);
      }
    }
  }
  if(err > tol)
  {
    printf("Maximum number of iterations reached.\n");
  }
  return g;
  
}

//' Matrix shrinkage by g-modeling method
//' 
//' The function use shrinkage method to estimate covariance matrix with given multinormal data. 
//' This method applies empirical method and approximate the optimal separable decision rule by
//' EM algorithm.
//' 
//' @param X          \eqn{n x p} matrix of data generated from p-variate multivariate normal distribution
//' @param d          interger. It controls the number of support points. Default value is 10.
//' @param centered   bool. Whether X has been centered. If center=false, it is known that the multivariate
//'                   distribution has mean of zeros. Otherwise means are unknown. Default value is false.
//' @param maxit      maximum number of allowable iterations. Default value is 200.
//' @param tol        error tolerance for convergence of EM algorithm. Default value is 1e-04.
//' @param verbose    bool. Whether print out error information during EM iteration. Default value is false.
//' @return 
//' 
//' - `est_cov`: \eqn{p x p} matrix of estiamted covariance matrix.
//' 
//' @examples
//' \donttest{
//' ## set parameters
//' p = 10
//' d = 20
//' n = 100
//' ## generate multivariate normal data X, with mean of zeros and identity matrix as covariance matrix.
//' X = matrix(rnorm(n*p,0,1), n, p)
//' ## estimate covariance with given data X
//' cov = msg(X=X, d=d)
//' }
//' 
//' @useDynLib cole
//' @export
// [[Rcpp::export]]
arma::mat msg(arma::mat X, int d = 10, bool centered = false, int maxit = 200, double tol = 0.0001, bool verbose=false)
{
  int n = X.n_rows;
  int p = X.n_cols;
  
  if(centered) n=n-1; // if the data matrix is centered, meaning means are unknown instead of zeros
  
  arma::mat est_cov(p,p);
  arma::mat S = trans(X) * X;
  
  // make grids
  arma::vec diag_entries(p);
  diag_entries = S.diag();
  arma::vec sds = sqrt(diag_entries / n);
  // build up sequences for standard deviations and correlations
  arma::vec grid_sd = linspace(sds.min(), sds.max(), d);
  arma::vec grid_rho = linspace(-1,1,d);
  grid_rho = grid_rho.subvec(1,d-2); //take out -1 and 1 from correlation sequence
  // build grid matrix for covariance. Three columns represent sd for row, sd for column, and correlation. 
  arma::mat grid(pow(d,2)*(d-2), 3); 
  grid.col(0) = rep_time(grid_sd, d*(d-2));
  grid.col(1) = rep_time(rep_each(grid_sd, d),d-2);
  grid.col(2) = rep_each(grid_rho, pow(d,2));
  grid = grid.rows(find(grid.col(0)<grid.col(1))); //remove half points according to the symmetry
  // build grid matrix for variance. 
  arma::mat grid_dg(d,3);
  grid_dg.col(0) = grid_sd;
  grid_dg.col(1) = grid_sd;
  grid_dg.col(2) = grid_sd.fill(1);
  // combine grid matrix for covariance and variance.
  grid = join_cols(grid, grid_dg);
  
  // position indices for lower triangular of S
  arma::mat cmbn0(p*p,2);
  arma::vec v = linspace(0,p-1,p);
  cmbn0.col(0) = rep_time(v,p);
  cmbn0.col(1) = rep_each(v,p);
  cmbn0 = cmbn0.rows(find(cmbn0.col(0) > cmbn0.col(1)));
  arma::umat cmbn = conv_to<arma::umat>::from(cmbn0);
  
  // build density matrix 
  mat W_all(cmbn.n_rows,4); //contains vectorised 2x2 matrices
  W_all.col(0) = S.elem(sub2ind(size(S),repmat(cmbn.col(0).t(),2,1)));
  W_all.col(1) = S.elem(sub2ind(size(S),cmbn.t()));
  W_all.col(2) = W_all.col(1);
  W_all.col(3) = S.elem(sub2ind(size(S),repmat(cmbn.col(1).t(),2,1)));
  arma::vec W_det = W_all.col(0) % W_all.col(3) - pow(W_all.col(1),2); //determinants of 2x2 matrices
  
  arma::vec log_det_w = log(W_det);
  arma::vec log_I1 = (n-3) / 2.0 * log_det_w;
  double log_I4 = n * log(2.0) + log(M_PI) / 2.0 + log(::tgamma(n/2.0)) + log(::tgamma((n-1)/2.0));
  
  // declare variables used in constructing density matrix
  arma::vec x(3);
  double prod1, log_I3;
  arma::vec sigma_v(4);
  arma::mat sigma(2,2);
  arma::mat sigma_inv(2,2);
  arma::vec log_I2(W_all.n_rows);
  arma::vec off(p*(p-1)/2), on(p);
  arma::vec tmp(p);
  arma::mat density_mat = mat(cmbn.n_rows+p,grid.n_rows); 
  for(unsigned int i_grid=0;i_grid<grid.n_rows;i_grid++)
  {
    x = grid.row(i_grid).t();
    if(x[2]<1)
    {
      prod1 = x[0]*x[1]*x[2];
      sigma = {{pow(x[0],2), prod1}, {prod1, pow(x[1],2)}};
      sigma_inv = sigma.i();
      
      log_I2 = W_all * vectorise(sigma_inv) / 2.0;
      log_I3 = n / 2.0 * log(det(sigma));
      off = exp(log_I1 - log_I2 - log_I3 - log_I4); //wishart density for covariances
      on = zeros(p);
      density_mat.col(i_grid) = join_cols(off,on);
    }
    else
    {
      off = zeros(p*(p-1)/2);
      tmp = S.diag() / pow(x[0],2);
      Rcpp::NumericVector on;
      on = dchisq(as<NumericVector>(wrap(tmp)),n) * n / pow(x[0], 2); //chi-square density for variances
      density_mat.col(i_grid) = join_cols(off,as<arma::vec>(wrap(on)));
    }
    
  }
  
  arma::uvec ind_inf = find(density_mat == R_PosInf);
  arma::vec new_inf = zeros(ind_inf.size());
  double c_max = density_mat.elem(find_finite(density_mat)).max();
  density_mat.elem(ind_inf) = new_inf.fill(c_max);
  // density_mat.elem( find_nonfinite(density_mat) ).fill(exp(-300));
  // determine the weights of grid points using EM algorithm
  arma::vec g = npmle_c(density_mat,maxit,tol,verbose); 
  
  arma::vec denom = density_mat * g;
  arma::vec tmp_cov(denom.size());
  arma::vec prods = prod(grid,1);
  tmp_cov = density_mat * (g % prods) / denom;
  arma::uvec vecind = sub2ind( size(est_cov), cmbn.t() );
  int p_cut = p*(p-1)/2;
  est_cov.elem(vecind) = tmp_cov.subvec(0,p_cut-1);
  est_cov.diag() = tmp_cov.subvec(p_cut, tmp_cov.size()-1);
  est_cov = symmatl(est_cov);
  
  return est_cov;
}
