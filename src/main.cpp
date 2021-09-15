
#include "main.h"
#include "misc_v11.h"
#include "probability_v11.h"

using namespace std;

Rcpp::List dummy1_cpp(Rcpp::List args) {
  
  vector<vector<double>> sigma = rcpp_to_matrix_double(args["sigma"]);
  vector<vector<double>> psi = rcpp_to_matrix_double(args["psi"]);
  double nu = rcpp_to_double(args["nu"]);
  
  double d1 = dinvwish2(sigma, psi, nu);
  print(d1);
  
  return Rcpp::List::create(Rcpp::Named("foo") = -9);
}

