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


//' Compound decision method for multivariate linear models
//'
//' The function uses EM algorithm to solve multivariate linear regression problems 
//' \deqn{Y = XB + \epsilon}
//' both outcome \eqn{Y} and feature \eqn{X} are multi-dimensional. Users can set distinct residual estimations for different outcomes or set identical estimation for more robust results. Details can be found in \href{https://github.com/sdzhao/cole}{our paper}.
//'
//' @param y \eqn{n x q} matrix of outcomes for training.
//' @param x \eqn{n x p} matrix of features for training. 
//' @param S \eqn{d x L} matrix of support points. If \eqn{L = p + 1}, then the first p columns are \eqn{\beta}s and the last column is the corresponding residual error estimates. If \eqn{L = p}, then each column of S is a vector of \eqn{\beta}s and argument min_s2 is required. \eqn{d = q x g} where g is the number of groups of support points. Support points can be estimated by other methods that solve multivariate linear regression. Eg. LASSO from glmnet, CMR from camel.
//' @param tol error tolerance for convergence of EM algorithm. Default value of tol is 1e-6.
//' @param maxit maximum number of allowable iterations. Default value of maxit is 1e5.
//' @param min_s2 a positive number corresponds to minimal variance of estimated y, min_s2 is required when there are p columns in S.
//'
//' @return
//' \item{f}{vector with \eqn{g x q} elements that describes the mixture of \eqn{\beta}s.}
//' \item{A}{matrix with dimension \eqn{q x d}. A is an estimation of likelihood by EM algorithm and will be used in predicting.}
//' \item{bs}{matrix with dimension \eqn{d x (p + 1)}. bs is support points used in updating prior distribution f. bs is equivalent to S when \eqn{L = p + 1} and will be used in predicting.}
//'
//' @examples
//' \donttest{
//' ## generate data
//' p = 10
//' q = 5
//' n = 50
//' x = matrix(rnorm(n*p,0,10), n, p)
//' beta = matrix(rnorm(p*q,0,10), q, p)
//' e = matrix(rnorm(n*q,0,0.1),n,q)
//' y = x %*% t(beta) + e
//' s2 = matrix(rep(0.1,q), q, 1)
//' ## initialize parameters for EM algorithm 
//' x_test = matrix(rnorm(n*p,0,1), n, p)
//' ## set minimal variance estimation min_s2 = 0.1
//' output1 = comte(y=y, x=x, S=beta, min_s2=0.1)
//' ## use distinct variance from multivariate linear regression models
//'output2 = comte(y=y, x=x, S=cbind(beta,s2))
//' }
//'
//' @useDynLib cole
//' @export
// [[Rcpp::export]]
List comte(arma::mat y, arma::mat x, arma::mat S, double tol = 0.000001 , int maxit = 100000, Nullable<NumericVector> min_s2 = R_NilValue){

  //initialization:
  int mp1 = S.n_rows; // g * q
  arma::vec f = arma::ones(mp1) / mp1;
  int p = x.n_cols;
  int q = y.n_cols;

  arma::mat B(mp1, p);
  arma::mat ss(mp1, 1);
  int p1 = S.n_cols;
  if(p1 == p){
    if(min_s2.isNull()){
      std::cout << "user defined minimum sigma squared is required" << std::endl;
    }else{
      B = S;
      NumericVector min_s2d(min_s2);
      ss.fill(min_s2d(0));
    }
  }else if(p1 == p + 1){
    B = S.submat(0,0,mp1-1,p-1);
    ss = S.submat(0,p,mp1-1,p);
  }else{
    std::cout << "wrong dimension" << std::endl;
  }

  //indices:
  arma::mat bs( mp1, p+1);
  bs = arma::join_rows(B,ss);

  //compute likelihood
  arma::mat ym = arma::repmat(y,1,mp1);
  arma::mat xbm = arma::kron(x * B.t(), arma::ones(1,q));
  arma::mat A(q, mp1);

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

  if(p1 == p + 1){
  	arma::mat ssm = arma::kron(ss.t(),arma::ones(1,q));
  	arma::mat Avec = exp( (-arma::sum(pow((ym - xbm),2), 0))/(2*ssm) ) + pow(0.1,300);
    A = arma::reshape(Avec, q, mp1);

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
      err = std::abs(diff(0,0)) / std::abs((oldll)(0,0));
      if(err <= tol){
        break;
      }
  	}
  }else{
  	arma::mat Avec = arma::sum(2*ym%xbm - xbm%xbm,0)/(2*ss(0,0));
  	float maxA = Avec.max();
  	Avec = exp(Avec / maxA);
  	A = arma::reshape(Avec, q, mp1);

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
  }

  return List::create(_["f"] = f, _["A"] = A, _["bs"] = bs);
}


//' Predicting function for compound decision method 
//'
//' The function takes returned object from comte and testing data as input. Return value is the prediction for corresponding output.
//'
//' @param comte_obj returned object from comte function that contains important information for predicting.
//' @param newx \eqn{m x p} matrix corresponds to testing data.
//'
//' @return
//' \item{esty}{\eqn{m x q} matrix of estimated y based on newx. }
//'
//' @examples
//' \donttest{
//' ## generate data
//' p = 10
//' q = 5
//' n = 50
//' x = matrix(rnorm(n*p,0,10), n, p)
//' beta = matrix(rnorm(p*q,0,10), q, p)
//' e = matrix(rnorm(n*q,0,0.1),n,q)
//' y = x %*% t(beta) + e
//' s2 = matrix(rep(0.1,q), q, 1)
//' ## initialize parameters for EM algorithm 
//' x_test = matrix(rnorm(n*p,0,1), n, p)
//' ## set minimal variance estimation min_s2 = 0.1
//' output1 = comte(y=y, x=x, S=beta, min_s2=0.1)
//' esty1 = predict_comte(output1, x_test)
//' ## use distinct variance from multivariate linear regression models
//'	output2 = comte(y=y, x=x, S=cbind(beta,s2))
//' esty2 = predict_comte(output2, x_test)
//' }
//'
//' @useDynLib cole
//' @export
// [[Rcpp::export]]
arma::mat predict_comte(List comte_obj, arma::mat newx){
  
  arma::vec f = comte_obj["f"];
  arma::mat A = comte_obj["A"];
  arma::mat bs = comte_obj["bs"];
  
  int q = A.n_rows;
  int mp1 = bs.n_rows;
  int n_new = newx.n_rows;
  int p = newx.n_cols;

  arma::mat esty(n_new,q);
  arma::mat fmat = arma::repmat(f,1,q);
  arma::mat numeritor2 = (newx * bs.submat(0,0,mp1-1,p-1).t()) * (A%fmat.t()).t() ;
  arma::mat denomiator2 = arma::repmat(1/(A*f),1, n_new);
  esty = numeritor2 % denomiator2.t();

  return esty;

}





