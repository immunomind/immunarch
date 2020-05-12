#include <Rcpp.h>
using namespace Rcpp;

// // [[Rcpp::export]]
// List rcpp_hello_world() {
//
//   CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
//   NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
//   List z            = List::create( x, y ) ;
//
//   return z ;
// }


// [[Rcpp::export]]
NumericVector fill_vec(NumericVector read_vec, NumericVector read_indices) {
  int dummy = 0;
  for (int i = 0; i < read_vec.size(); i++) {
    for (int j = dummy; j < (read_vec[i] + dummy); j++) {
      read_indices[j] = i;
    }
    dummy = dummy + read_vec[i];
  }
  return read_indices;
}

// [[Rcpp::export]]
NumericVector fill_reads(NumericVector new_reads, NumericVector new_counts) {
  for (int i = 0; i < new_counts.size(); i++) {
    new_reads[new_counts[i]] = new_reads[new_counts[i]] + 1;
  }
  return new_reads;
}
