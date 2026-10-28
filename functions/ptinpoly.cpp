#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' @param points A \code{rowvec} with x,y  coordinate structure.
//' @param bp     A \code{matrix} containing the boundary points of the polygon. 
//' @return A \code{bool} indicating whether the point is in the polygon (TRUE) or not (FALSE)
// [[Rcpp::export]]
bool pnpoly(const arma::rowvec& point, const arma::mat& bp) {
    // Implementation of the ray-casting algorithm is based on
    // 
    unsigned int i, j;
    
    double x = point(0), y = point(1);
    
    bool inside = false;
    for (i = 0, j = bp.n_rows - 1; i < bp.n_rows; j = i++) {
      double xi = bp(i,0), yi = bp(i,1);
      double xj = bp(j,0), yj = bp(j,1);
      
      // See if point is inside polygon
      inside ^= (((yi >= y) != (yj >= y)) && (x <= (xj - xi) * (y - yi) / (yj - yi) + xi));
    }
    
    // Is the cat alive or dead?
    return inside;
}