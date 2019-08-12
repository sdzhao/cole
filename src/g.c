/* ***************************************************************
g.c
functions for single gaussian sequence problem
*************************************************************** */
#include <R.h>
#include <Rinternals.h>
#include <math.h>

// N(0,1) density
double phi(double x){
  return(1 / sqrt(2 * M_PI) * exp(-0.5 * (x * x)));
}

// helper function to calculate weighted sums
double * g_wsum(double x1, double s1, int n,
		double *t1, double rho){
  static double ret[4];
  
  // initialize
  ret[0] = 0; ret[1] = 0; ret[2] = 0;
  
  // calculate
  double tmp;
  for(int j = 0; j < n; j++){
    tmp = 1 / s1 * phi((x1 - t1[j]) / s1);
    ret[0] += tmp;
    ret[1] += (t1[j] - x1) * tmp;
    ret[2] += (t1[j] - x1) * (t1[j] - x1) * tmp;
  }
  
  return(ret);
}

// regularized nonparametric separable rule
SEXP g_rule(SEXP R_x1, SEXP R_s1,
	    SEXP R_t1, SEXP R_rho){
  int n = length(R_x1);
  double *x1, *t1;
  x1 = REAL(R_x1);
  t1 = REAL(R_t1);
  double s1, rho;
  s1 = REAL(R_s1)[0];
  rho = REAL(R_rho)[0];
  
  SEXP R_rule = PROTECT(allocVector(REALSXP, n));
  double *rule = REAL(R_rule);

  double * wsum;
  for(int i = 0; i < n; i++){
    wsum = g_wsum(x1[i], s1, n, t1, rho);
    rule[i] = x1[i] + wsum[1] / (rho + wsum[0]);
  }
  
  UNPROTECT(1);
  return(R_rule);
}

// SURE
SEXP g_sure(SEXP R_x1, SEXP R_s1,
	    SEXP R_t1, SEXP R_rho){
  int n = length(R_x1);
  double *x1, *t1;
  x1 = REAL(R_x1);
  t1 = REAL(R_t1);
  double s1, rho;
  s1 = REAL(R_s1)[0];
  rho = REAL(R_rho)[0];
  
  SEXP R_sure = PROTECT(allocVector(REALSXP, 1));
  double *sure = REAL(R_sure);
  // initialize
  sure[0] = 0;
  
  double * wsum;
  for(int i = 0; i < n; i++){
    wsum = g_wsum(x1[i], s1, n, t1, rho);
    sure[0] += 2 * wsum[2] / (rho + wsum[0]) -
      2 * s1 * s1 * wsum[0] / (rho + wsum[0]) -
      wsum[1] * wsum[1] / (rho + wsum[0]) / (rho + wsum[0]);
  }
  
  sure[0] = sure[0] / n + s1 * s1;

  UNPROTECT(1);
  return(R_sure);
}

// minimize sure
SEXP g_min_sure(SEXP R_x1, SEXP R_s1,
		SEXP R_rho,
		SEXP R_K, SEXP R_C,
		SEXP R_tol, SEXP R_maxit, SEXP R_maxbreak, SEXP R_verbose){
  int n = length(R_x1);
  int i,j; // i indexes x's, j's index t's
  
  double *x1 = REAL(R_x1);
  double s1, rho;
  s1 = REAL(R_s1)[0];
  rho = REAL(R_rho)[0];
  int K = INTEGER(R_K)[0];
  double C = REAL(R_C)[0];
  double tol = REAL(R_tol)[0];
  int maxit = INTEGER(R_maxit)[0];
  int maxbreak = INTEGER(R_maxbreak)[0];
  int verbose = INTEGER(R_verbose)[0];
  
  SEXP R_t = PROTECT(allocVector(REALSXP, n));
  double *t = REAL(R_t);
  // initialize ts
  for(j = 0; j < n; j++){
      t[j] = x1[j];
  }
  
  // precalculate quantities in the SURE formula
  double wsum0[n], wsum1[n], wsum2[n], * wsum;
  for(i = 0; i < n; i++){
    wsum = g_wsum(x1[i], s1, n, x1, rho);
    wsum0[i] = wsum[0];
    wsum1[i] = wsum[1];
    wsum2[i] = wsum[2];
  }
  
  // initial sure value
  double sure_old = 0;
  for(i = 0; i < n; i++){
    sure_old += 2 * wsum2[i] / (rho + wsum0[i]) -
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
  double t1_new, inc; // also try t[j]s
  double tmp, wsum0_new[n], wsum1_new[n], wsum2_new[n];
  double sure_min, sure_new; // sure values when trying t[j]s
  int optk; // optimal index k of the new t1/t2[j]
  
  // don't break at first convergence of SURE, in case going back
  // to change beginning coordinates drops SURE by a lot
  int breakcount = 0;
  sure_min = sure_old;
  for(iter = 0; iter < maxit; iter++){
    
    // coordinate descent
    for(j = 0; j < n; j++){
      
      // lower and upper bounds
      l = x1[j] - C * s1;
      u = x1[j] + C * s1;
      
      // subtract current t1 and t2 terms
      for(i = 0; i < n; i++){
	tmp = 1 / s1 * phi((x1[i] - t[j]) / s1);
	wsum0[i] -= tmp;
	wsum1[i] -= (t[j] - x1[i]) * tmp;
	wsum2[i] -= (t[j] - x1[i]) * (t[j] - x1[i]) * tmp;
      }
      
      // determine optimal value of t[j]
      // consider original value and K equally spaced values in [l, u]
      // if k = 0 then try original t[j] value
      inc = (u - l) / (K - 1);
      optk = 0;
      for(k = 0; k < (K + 1); k++){
	
	// new value of t1 to try
	if(k < 1){
	  t1_new = t[j];
	} else {
	  t1_new = l + (k - 1) * inc;
	}
	  
	// update values of precalculated quantities
	for(i = 0; i < n; i++){
	  tmp = 1 / s1 * phi((x1[i] - t1_new) / s1);
	  wsum0_new[i] = wsum0[i] + tmp;
	  wsum1_new[i] = wsum1[i] + (t1_new - x1[i]) * tmp;
	  wsum2_new[i] = wsum2[i] + (t1_new - x1[i]) * (t1_new - x1[i]) * tmp;
	}
	  
	// calculate sure
	sure_new = 0;
	for(i = 0; i < n; i++){
	  sure_new += 2 * wsum2_new[i] / (rho + wsum0_new[i]) -
	    2 * s1 * s1 * wsum0_new[i] / (rho + wsum0_new[i]) -
	    wsum1_new[i] * wsum1_new[i] / (rho + wsum0_new[i]) / (rho + wsum0_new[i]);
	}
	sure_new = sure_new / n + s1 * s1;

	// update optimal indices
	if(sure_new <= sure_min){
	  optk = k;
	  sure_min = sure_new;
	}
	
      } // end loop over k, for determining optimal value of t[j]
      
      // update t[j] and wsum's
      if(optk > 0){
	t[j] = l + (optk - 1) * inc;
      }
      for(i = 0; i < n; i++){
	tmp =  1 / s1 * phi((x1[i] - t[j]) / s1);
	wsum0[i] += tmp;
	wsum1[i] += (t[j] - x1[i]) * tmp;
	wsum2[i] += (t[j] - x1[i]) * (t[j] - x1[i]) * tmp;
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
