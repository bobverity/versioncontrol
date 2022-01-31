
#include "main.h"
#include "misc_v12.h"
#include "probability_v14.h"
#include "Sampler_v4.h"

using namespace std;

Rcpp::List dummy1_cpp(Rcpp::List args) {
  
  int n = rcpp_to_int(args["n"]);
  double alpha = rcpp_to_double(args["alpha"]);
  vector<double> z = rdirichlet1(alpha, n);
  
  return Rcpp::List::create(Rcpp::Named("z") = z);
}

