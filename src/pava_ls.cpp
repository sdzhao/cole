#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

//' Least-Squares Pooled Adjacent Violators Algorithm
//'
//' An implementation of the pooled adjacent violators algorithm for
//' minimizing the least squares loss subject to the monotonicity
//' (non-decreasing) constraint.
//'
//' It is assumed that `y` is ordered corresponding to the correct indices for
//' a monotone non-decreasing fit. If a monotone non-increasing fit is desired
//' use `base::rev()` as in the examples.
//'
//' The basic idea is to start by partitioning the data into singletons (zero
//' least squares loss), then merge neighboring sets if the monotone
//' non-decreasing constraint is violated. The value a set takes is the average
//' of its elements, this is implemented by tracking set cardinality and
//' the current mean so that the merger of two sets is just a weighted sum.
//'
//' This algorithm has time complexity O(n). However, this implementation may
//' have worse complexity because of the way we are removing elements in a
//' vector during the merge.
//'
//' @param y numeric vector of observations corresponding to a sorted,
//' increasing predictor
//' @return a  numeric vector of fitted values corresponding to y
//' @examples
//' # non-decreasing example
//' set.seed(1)
//' y <- atan(seq(from = -5, to = 5, length.out = 51)) + rnorm(51, sd = 0.2)
//' y.hat <- pava_ls(y) # not exported
//' plot(y)
//' points(y.hat, col = "red")
//'
//' #non-increasing example
//' set.seed(1)
//' y <- -atan(seq(from = -5, to = 5, length.out = 51)) + rnorm(51, sd = 0.2)
//' y.hat <- rev(pava_ls(rev(y)))
//' plot(y)
//' points(y.hat, col = "red")
//'
//' @useDynLib cole
//' @importFrom Rcpp evalCpp
//' @export
// [[Rcpp::export]]
NumericVector pava_ls(NumericVector y) {
  // initialization
  int n = y.length();
  int m = 0;
  std::vector<int> blocksize(n, 1);
  std::vector<double> values(y.begin(), y.end());

  // fit
  int i = 0;
  while (i < n - m - 1) {
    // check merge condition
    if (values[i] >= values[i+1]) {
      // merge
      values[i] = (blocksize[i]*values[i] + blocksize[i+1]*values[i+1]) / (blocksize[i]+blocksize[i+1]);
      blocksize[i] = blocksize[i] + blocksize[i+1];

      // drop redundant bins
      values.erase(values.begin() + i + 1);
      blocksize.erase(blocksize.begin() + i + 1);
      m = m + 1;

      // step back
      if (i > 0) {
        i = i - 1;
      }
    }
    else {
      i = i + 1;
    }
  }

  // expand blocks and values to correspond to y
  NumericVector out(n);
  int idx = 0;
  for (int i = 0; i < blocksize.size(); i++) {
    for (int j = 0; j < blocksize[i]; j ++) {
      out[idx] = values[i];
      idx++;
    }
  }

  return out;
}
