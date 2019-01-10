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


// [[Rcpp::export]]
List comte(arma::mat y, arma::mat x, arma::mat S, double tol , int maxit, int p, int q, int n, arma::mat newx, Nullable<NumericVector> min_s2 = R_NilValue){

  //initialization:
  int mp1 = S.n_rows; // g * q
  arma::vec f = arma::ones(mp1) / mp1;

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
    err = std::abs(diff(0,0)) / std::abs((oldll)(0,0));
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
