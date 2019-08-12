/* ***************************************************************
gg.c
functions for gaussian sequence problem with gaussian side information
*************************************************************** */
#include <R.h>
#include <Rinternals.h>
#include <math.h>

// bivariate N(mu, S) density
double phi2(double x1, double u1, double s1,
	    double x2, double u2, double s2,
	    double r){
  double z = (double) (x1 - u1) * (x1 - u1) / s1 / s1 -
    2 * r * (x1 - u1) * (x2 - u2) / s1 / s2 +
    (x2 - u2) * (x2 - u2) / s2 / s2;
  return(1 / (2 * M_PI * s1 * s2 * sqrt(1 - r * r)) *
	 exp(-0.5 * z / (1 - r * r)));
}

// helper function to calculate weighted sums
double * gg_wsum(double x1, double s1, double x2, double s2,
		 double r,
		 int n,
		 double *t, double rho){
  static double ret[3];
  
  // initialize
  ret[0] = 0; ret[1] = 0; ret[2] = 0;
  
  // calculate
  double tmp;
  for(int j = 0; j < n; j++){
    tmp = phi2(x1, t[j], s1, x2, t[j + n], s2, r);
    ret[0] += tmp;
    ret[1] += (t[j] - x1) * tmp;
    ret[2] += (t[j] - x1) * (t[j] - x1) * tmp;
  }
  
  return(ret);
}

// regularized nonparametric separable rule
SEXP gg_rule(SEXP R_x1, SEXP R_s1, SEXP R_x2, SEXP R_s2,
	     SEXP R_r,
	     SEXP R_t, SEXP R_rho){
  int n = length(R_x1);
  double *x1, *x2, *t;
  x1 = REAL(R_x1);
  x2 = REAL(R_x2);
  t = REAL(R_t);
  double s1, s2, r, rho;
  s1 = REAL(R_s1)[0];
  s2 = REAL(R_s2)[0];
  r = REAL(R_r)[0];
  rho = REAL(R_rho)[0];
  
  SEXP R_rule = PROTECT(allocVector(REALSXP, n));
  double *rule = REAL(R_rule);

  double * wsum;
  for(int i = 0; i < n; i++){
    wsum = gg_wsum(x1[i], s1, x2[i], s2, r, n, t, rho);
    rule[i] = x1[i] + wsum[1] / (rho + wsum[0]);
  }
  
  UNPROTECT(1);
  return(R_rule);
}

// SURE
double gg_sure(double *x1, double s1, double *x2, double s2,
	       int n, double *t, double rho) {

  double sure = 0;
  
  double * wsum;
  for(int i = 0; i < n; i++){
    wsum = gg_wsum(x1[i], s1, x2[i], s2, 0, n, t, rho);
    sure += 2 * wsum[2] / (rho + wsum[0]) +
      2 * s1 * s1 * wsum[0] / (rho + wsum[0]) -
      wsum[1] * wsum[1] / (rho + wsum[0]) / (rho + wsum[0]);
  }
  
  sure = sure / n + s1 * s1;

  return(sure);
  
}

// fast algorithm adds/subtracts potentially small numbers
// and can result in numeric instability, so requires
// rho > 0, e.g. rho = 1e-12
// minimize sure one coordinate at a time
SEXP gg_min_sure_fast(SEXP R_x1, SEXP R_s1, SEXP R_x2, SEXP R_s2,
		      SEXP R_rho,
		      SEXP R_K, SEXP R_C,
		      SEXP R_tol, SEXP R_maxit, SEXP R_maxbreak,
		      SEXP R_verbose){
  int n = length(R_x1);
  int i,j; // i indexes x's, j's index t's
  
  double *x1, *x2;
  x1 = REAL(R_x1);
  x2 = REAL(R_x2);
  double s1, s2, rho;
  s1 = REAL(R_s1)[0];
  s2 = REAL(R_s2)[0];
  rho = REAL(R_rho)[0];
  int K = INTEGER(R_K)[0];
  double C = REAL(R_C)[0];
  double tol = REAL(R_tol)[0];
  int maxit = INTEGER(R_maxit)[0];
  int maxbreak = INTEGER(R_maxbreak)[0];
  int verbose = INTEGER(R_verbose)[0];
  
  // first n components are t1, second n are t2
  SEXP R_t = PROTECT(allocVector(REALSXP, 2 * n));
  double *t = REAL(R_t);
  // initialize ts
  for(j = 0; j < (2 * n); j++){
    if(j < n){
      t[j] = x1[j];
    } else {
      t[j] = x2[j - n];
    }
  }
  
  // precalculate quantities in the SURE formula
  double wsum0[n], wsum1[n], wsum2[n], * wsum;
  for(i = 0; i < n; i++){
    wsum = gg_wsum(x1[i], s1, x2[i], s2, 0, n, t, rho);
    wsum0[i] = wsum[0];
    wsum1[i] = wsum[1];
    wsum2[i] = wsum[2];
  }

  // initial sure value
  double sure_old = 0;
  for(i = 0; i < n; i++){
    sure_old += 2 * wsum2[i] / (rho + wsum0[i]) +
      2 * s1 * s1 * wsum0[i] / (rho + wsum0[i]) -
      wsum1[i] * wsum1[i] / (rho + wsum0[i]) / (rho + wsum0[i]);
  }
  sure_old = sure_old / n + s1 * s1;
  if(verbose == 1){
    Rprintf("initial SURE: %f\n", sure_old);
  }
  
  // find minimizing t
  int iter, k;
  double l, u; // lower and upper bounds
  double t1_new, t2_new, inc; // also try t[j]s
  double tmp, wsum0_new[n], wsum1_new[n], wsum2_new[n];
  double sure_min, sure_new; // sure values when trying t[j]s
  int optk; // optimal index k of the new t1/t2[j]
  
  // don't break at first convergence of SURE, in case going back
  // to change beginning coordinates drops SURE by a lot
  int breakcount = 0;
  sure_min = sure_old;
  for(iter = 0; iter < maxit; iter++){
    
    // coordinate descent (coordinate NOT by pair)
    for(j = 0; j < (2 * n); j++){
      
      // lower and upper bounds
      if(j < n){ // updating t1 coordinates
	t1_new = t[j]; // current t1
	t2_new = t[j + n]; // current t2
	// sd = s1
	l = x1[j] - C * s1;
	u = x1[j] + C * s1;
      } else { // updating t2 coordinates
	t1_new = t[j - n]; // current t1
	t2_new = t[j]; // current t2
	// sd = s2
	l = x2[j - n] - C * s2;
	u = x2[j - n] + C * s2;
      }
      
      // subtract current t1 and t2 terms
      for(i = 0; i < n; i++){
	tmp = phi2(x1[i], t1_new, s1, x2[i], t2_new, s2, 0);
	wsum0[i] -= tmp;
	wsum1[i] -= (t1_new - x1[i]) * tmp;
	wsum2[i] -= (t1_new - x1[i]) * (t1_new - x1[i]) * tmp;
      }
      
      // determine optimal value of t[j]
      // consider original value and K equally spaced values in [l, u]
      // if k = 0 then try original t[j] value
      inc = (u - l) / (K - 1);
      optk = 0;
      for(k = 0; k < (K + 1); k++){

	if(j < n){
	  // new value of t1 to try
	  if(k < 1){
	    t1_new = t[j];
	  } else {
	    t1_new = l + (k - 1) * inc;
	  }
	  // keep value of t2_new the same
	} else {
	  // keep value of t1_new the same
	  // new value of t2 to try
	  if(k < 1){
	    t2_new = t[j];
	  } else {
	    t2_new = l + (k - 1) * inc;
	  }
	}
	  
	// update values of precalculated quantities
	for(i = 0; i < n; i++){
	  tmp = phi2(x1[i], t1_new, s1, x2[i], t2_new, s2, 0);
	  wsum0_new[i] = wsum0[i] + tmp;
	  wsum1_new[i] = wsum1[i] + (t1_new - x1[i]) * tmp;
	  wsum2_new[i] = wsum2[i] + (t1_new - x1[i]) * (t1_new - x1[i]) * tmp;
	}
	
	// calculate sure
	sure_new = 0;
	for(i = 0; i < n; i++){
	  sure_new += 2 * wsum2_new[i] / (rho + wsum0_new[i]) +
	    2 * s1 * s1 * wsum0_new[i] / (rho + wsum0_new[i]) -
	    wsum1_new[i] * wsum1_new[i] / (rho + wsum0_new[i]) / (rho + wsum0_new[i]);
	}
	sure_new = sure_new / n + s1 * s1;
	
	if(sure_new <= sure_min){
	  optk = k;
	  sure_min = sure_new;
	}
	  
      } // end loop over k, for determining optimal value of t[j]
      
      // update t[j] and wsum's
      if(optk > 0){
	t[j] = l + (optk - 1) * inc;
      }
      if(j < n){
	t1_new = t[j];
	t2_new = t[j + n];
      } else {
	t1_new = t[j - n];
	t2_new = t[j];
      }
      for(i = 0; i < n; i++){
	tmp = phi2(x1[i], t1_new, s1, x2[i], t2_new, s2, 0);
	wsum0[i] += tmp;
	wsum1[i] += (t1_new - x1[i]) * tmp;
	wsum2[i] += (t1_new - x1[i]) * (t1_new - x1[i]) * tmp;
      }

    } // end coordinate descent
    
    if(verbose == 1){
      Rprintf("iter:%d,sure:%f\n", iter + 1, sure_min);
    }
    // check convergence
    // sure_new from the last coordinate is guaranteed to be
    // i)  lowest among all coordinates
    // ii) lower than or equal to sure_old
    if((sure_old - sure_min) < tol){
      breakcount++;
    }
    if(breakcount >= maxbreak){
      break;
    } else {
      sure_old = sure_min;
    }
    
    // if never break, warn that maxit was reached
    if(iter == (maxit - 1)){
      Rprintf("maxit reached\n");
    }
    
  } // end minimization of sure
  
  UNPROTECT(1);
  return(R_t);
}

// minimize sure one coordinate at a time
SEXP gg_min_sure(SEXP R_x1, SEXP R_s1, SEXP R_x2, SEXP R_s2,
		 SEXP R_rho,
		 SEXP R_K, SEXP R_C,
		 SEXP R_tol, SEXP R_maxit, SEXP R_maxbreak, SEXP R_verbose){
  int n = length(R_x1);
  int j; // j's index t's
  
  double *x1, *x2;
  x1 = REAL(R_x1);
  x2 = REAL(R_x2);
  double s1, s2, rho;
  s1 = REAL(R_s1)[0];
  s2 = REAL(R_s2)[0];
  rho = REAL(R_rho)[0];
  int K = INTEGER(R_K)[0];
  double C = REAL(R_C)[0];
  double tol = REAL(R_tol)[0];
  int maxit = INTEGER(R_maxit)[0];
  int maxbreak = INTEGER(R_maxbreak)[0];
  int verbose = INTEGER(R_verbose)[0];
  
  // first n components are t1, second n are t2
  SEXP R_t = PROTECT(allocVector(REALSXP, 2 * n));
  double *t = REAL(R_t);
  // initialize ts
  for(j = 0; j < (2 * n); j++){
    if(j < n){
      t[j] = x1[j];
    } else {
      t[j] = x2[j - n];
    }
  }

  // initial sure value
  double sure_old = gg_sure(x1, s1, x2, s2, n, t, rho);
  if(verbose == 1){
    Rprintf("initial SURE: %f\n", sure_old);
  }
  
  // find minimizing t
  int iter, k;
  double l, u; // lower and upper bounds
  double t_old, inc; // old value of t[j] if new one doesn't lower sure
  double sure_min, sure_new; // sure values when trying t[j]s
  
  // don't break at first convergence of SURE
  // need an extra cycle of no improvement in SURE to converge
  int breakcount = 0;
  sure_min = sure_old;
  for(iter = 0; iter < maxit; iter++){
    
    // coordinate descent (coordinate NOT by pair)
    for(j = 0; j < (2 * n); j++){
      
      // store old t[j] and calculate lower and upper bounds
      t_old = t[j];
      if(j < n){ // t1
	l = x1[j] - C * s1;
	u = x1[j] + C * s1;
      } else { // t2
	l = x2[j - n] - C * s2;
	u = x2[j - n] + C * s2;
      }
      
      // determine optimal value of t[j]
      // consider K equally spaced values in [l, u]
      inc = (u - l) / (K - 1);
      for(k = 0; k < K; k++){
	
	// try new value of t[j]
	t[j] = l + k * inc;

	// calculate sure
	sure_new = gg_sure(x1, s1, x2, s2, n, t, rho);
	
	// update t[j] only if sure is reduced
	if (sure_new <= sure_min) {
	  sure_min = sure_new;
	  t_old = t[j];
	} else {
	  t[j] = t_old;
	}
	  
      } // end loop over k, for determining optimal value of t[j]
      
    } // end coordinate descent
    
    if(verbose == 1){
      Rprintf("iter:%d,sure:%f\n", iter + 1, sure_min);
    }
    // check convergence
    // sure_new from the last coordinate is guaranteed to be
    // i)  lowest among all coordinates
    // ii) lower than or equal to sure_old
    if((sure_old - sure_min) < tol){
      breakcount++;
    }
    if(breakcount >= maxbreak){
      break;
    } else {
      sure_old = sure_min;
    }
    
    // if never break, warn that maxit was reached
    if(iter == (maxit - 1)){
      Rprintf("maxit reached\n");
    }
    
  } // end minimization of sure
  
  UNPROTECT(1);
  return(R_t);
}
