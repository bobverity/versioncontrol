
#pragma once

// comment out this definition to switch from Rcpp to C++
#define RCPP_ACTIVE

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

#include <iostream>
#include <vector>

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List dummy1_cpp(Rcpp::List args);
