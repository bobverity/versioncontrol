
#include "main.h"
#include "misc_v12.h"
#include "probability_v14.h"
#include "Sampler_v4.h"

using namespace std;

Rcpp::List dummy1_cpp(Rcpp::List args) {
  
  vector<double> x = rcpp_to_vector_double(args["x"]);
  vector<double> mu = rcpp_to_vector_double(args["mu"]);
  vector<vector<double>> sigma = rcpp_to_matrix_double(args["sigma"]);
  
  double z = dmnorm2(x, mu, sigma);
  print(z);
  
  return Rcpp::List::create(Rcpp::Named("foo") = -9);
}

