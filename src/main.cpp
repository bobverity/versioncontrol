
#include "main.h"
#include "misc_v13.h"
#include "probability_v17.h"
#include "Sampler_v4.h"

using namespace std;

Rcpp::List dummy1_cpp(Rcpp::List args) {
  
  double z = dbeta1(0.3, 2, 3);
  
  return Rcpp::List::create(Rcpp::Named("z") = z);
}

