// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

using namespace Rcpp;

//' Compound decision method for two-sieded models
//'
//' EM algorithm for predicting unseen data for two-sided linear model based on candidate parameters estimated by simple models
//'
//'
//' @param y outcomes for training with dimension n by q
//' @param x features for training with dimension n by p
//' @param ss condidates sigma square estimated by other mithods (lasso, ridge, camel)
//' @param B condidates beta estimated by other mithods (lasso, ridge, camel)
//' @param tol error tolerance for convergence
//' @param maxit maximum number of allowable iterations
//' @param p number of features
//' @param q number of outcomes
//' @param n number of observations for training
//' @param newx new observation for prediction
//'
//' @return
//' \item{f}{estimation of prior of beta}
//' \item{A}{likelihood of training data}
//' \item{bs}{candidates beta}
//' \item{esty}{estimation for based on newx}
//'
//' @examples
//' \donttest{
//' p = 10
//' q = 5
//' n = 50
//'
//' x = matrix(rnorm(n*p,0,10), n, p)
//' beta = matrix(rnorm(p*q,0,10), q, p)
//' e = matrix(rnorm(n*q,0,0.1),n,q)
//' y = x %*% t(beta) + e
//' s2 = matrix(rep(0.1,q), q, 1)
//' tol = 0.001
//' maxit = 1000
//' x_test = matrix(rnorm(n*p,0,1), n, p)
//'
//' output = comte(y, x, s2, beta, tol, maxit, p, q, n, x_test)
//'
//' }
//' @export
// [[Rcpp::export]]
List comte(arma::mat y, arma::mat x, arma::mat ss, arma::mat B, double tol , int maxit, int p, int q, int n, arma::mat newx){

  //initialization:
  int mp1 = B.n_rows; // g * q
  arma::vec f = arma::ones(mp1) / mp1;

  //indices:

  arma::mat bs( mp1, p+1);
  bs = arma::join_rows(B,ss);

  // //compute likelihood
  arma::mat ym = arma::repmat(y,1,mp1);
  arma::mat xbm = arma::kron(x * B.t(), arma::ones(1,q));
  arma::mat ssm = arma::kron(ss.t(),arma::ones(1,q));
  //arma::mat Avec = exp( (-arma::sum(pow((ym - xbm),2), 0))/(2*ssm))/pow(ssm, n/2) + pow(0.1,300);
  arma::mat Avec = exp( (-arma::sum(pow((ym - xbm),2), 0))/(2*ssm) ) + pow(0.1,300);
  arma::mat A = arma::reshape(Avec, q, mp1);

  //EM algorithm
  arma::mat ll;
  arma::mat oldll;
  arma::mat diff;
  arma::mat thres;

  arma::mat denomiator;
  arma::mat weight;

  double err = pow(10,10);
  arma::vec oldf;
  int tol_iter = 0;
  for(int i = 0; i < maxit; i++){
    tol_iter += 1;
    oldf = f;
    thres = 1/(A * oldf);
    for(int j = 0; j < q; j++){
      if( thres(j,0) > pow(10,300)){
        thres(j,0) = pow(10,300);
      }
    }


    f = A.t() * (thres) % oldf /q;
    ll = sum(log(A * f));
    oldll = sum(log(A * oldf));
    diff = oldll - ll;
    err = abs(diff(0,0)) / abs((oldll)(0,0));
    if(err <= tol){
      break;
    }
  }

  int n_new = newx.n_rows;
  arma::mat esty(n_new,q);
  arma::mat fmat = arma::repmat(f,1,q);
  arma::mat numeritor2 = (newx * bs.submat(0,0,mp1-1,p-1).t()) * (A%fmat.t()).t() ;
  arma::mat denomiator2 = arma::repmat(1/(A*f),1, n_new);
  esty = numeritor2 % denomiator2.t();

  return List::create(_["f"] = f, _["A"] = A, _["bs"] = bs, _["esty"] = esty, _["minA"] = A.min(),  _["maxA"] = A.max());
}

//' Compound decision method for two-sieded models with minimal sigma square
//'
//' EM algorithm for predicting unseen data for two-sided linear model based on candidate parameters estimated by simple models. Sigma square is taken to be the minimum for better accuracy in some situations
//'
//'
//' @param y outcomes for training with dimension n by q
//' @param x features for training with dimension n by p
//' @param ss condidates sigma square estimated by other mithods (lasso, ridge, camel)
//' @param B condidates beta estimated by other mithods (lasso, ridge, camel)
//' @param tol error tolerance for convergence
//' @param maxit maximum number of allowable iterations
//' @param p number of features
//' @param q number of outcomes
//' @param n number of observations for training
//' @param newx new observation for prediction
//'
//' @return
//' \item{f}{estimation of prior of beta}
//' \item{A}{likelihood of training data}
//' \item{bs}{candidates beta}
//' \item{esty}{estimation for based on newx}
//'
//' @examples
//' \donttest{
//' p = 10
//' q = 5
//' n = 50
//'
//' x = matrix(rnorm(n*p,0,10), n, p)
//' beta = matrix(rnorm(p*q,0,10), q, p)
//' e = matrix(rnorm(n*q,0,0.1),n,q)
//' y = x %*% t(beta) + e
//' s2 = matrix(rep(0.1,q), q, 1)
//' tol = 0.001
//' maxit = 1000
//' x_test = matrix(rnorm(n*p,0,1), n, p)
//'
//' output = comte_min(y, x, s2, beta, tol, maxit, p, q, n, x_test)
//'
//' }
//'
//' @export
// [[Rcpp::export]]
List comte_min(arma::mat y, arma::mat x, arma::mat ss, arma::mat B, double tol , int maxit, int p, int q, int n, arma::mat newx){

  //initialization:
  int mp1 = B.n_rows;
  arma::vec f = arma::ones(mp1) / mp1;

  //indices:

  arma::mat bs( mp1, p+1);
  bs = arma::join_rows(B,ss);

  // //compute likelihood
  arma::mat ym = arma::repmat(y,1,mp1);
  arma::mat xbm = arma::kron(x * B.t(), arma::ones(1,q));
  arma::mat Avec = arma::sum(2*ym%xbm - xbm%xbm,0)/(2*ss(0,0));
  float maxA = Avec.max();
  Avec = exp(Avec / maxA);
  arma::mat A = arma::reshape(Avec, q, mp1);

  //EM algorithm
  arma::mat ll;
  arma::mat oldll;
  arma::mat diff;
  arma::mat thres;

  arma::mat denomiator;
  arma::mat weight;

  double err = pow(10,10);
  arma::vec oldf;
  int tol_iter = 0;
  for(int i = 0; i < maxit; i++){
    tol_iter += 1;
    oldf = f;
    thres = 1/(A * oldf);
    f = A.t() * (thres) % oldf /q;
    ll = sum(log(A * f));
    oldll = sum(log(A * oldf));
    diff = oldll - ll;
    err = abs(diff(0,0)) / abs((oldll)(0,0));
    if(err <= tol){
      break;
    }
  }

  int n_new = newx.n_rows;
  arma::mat esty(n_new,q);
  arma::mat fmat = arma::repmat(f,1,q);
  arma::mat numeritor2 = (newx * bs.submat(0,0,mp1-1,p-1).t()) * (A%fmat.t()).t() ;
  arma::mat denomiator2 = arma::repmat(1/(A*f),1, n_new);
  esty = numeritor2 % denomiator2.t();

  return List::create(_["f"] = f, _["A"] = A, _["bs"] = bs, _["esty"] = esty);
}
