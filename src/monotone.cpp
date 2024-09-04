#include <stdlib.h>
#include <Rcpp.h>

void monotoneC( int *n, double* x, double* w )
// Function monotone(),
// performs simple linear ordered monotone regression
// Copyright (C) 2020 Frank M.T.A. Busing (e-mail: busing at fsw dot leidenuniv dot nl)
// This function is free software:
// you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// You should have received a copy of the GNU General Public License along with this function.
// If not, see <https://www.gnu.org/licenses/>.
{
  double* rx = &x[-1];
  double* rw = &w[-1];
  size_t* idx = ( size_t* ) calloc( *n + 1, sizeof( size_t ) );
  idx[0] = 0;
  idx[1] = 1;
  size_t b = 1;
  double xbm1 = rx[b];
  double wbm1 = rw[b];
  for ( size_t i = 2; i <= *n; i++ ) {
    b++;
    double xb = rx[i];
    double wb = rw[i];
    if ( xbm1 > xb ) {
      b--;
      double sb = wbm1 * xbm1 + wb * xb;
      wb += wbm1;
      xb = sb / wb;
      while ( i < *n && xb >= rx[i + 1] ) {
        i++;
        sb += rw[i] * rx[i];
        wb += rw[i];
        xb = sb / wb;
      }
      while ( b > 1 && rx[b - 1] > xb ) {
        b--;
        sb += rw[b] * rx[b];
        wb += rw[b];
        xb = sb / wb;
      }
    }
    rx[b] = xbm1 = xb;
    rw[b] = wbm1 = wb;
    idx[b] = i;
  }
  size_t from = *n;
  for ( size_t k = b; k > 0; k-- ) {
    const size_t to = idx[k - 1] + 1;
    const double xk = rx[k];
    for ( size_t i = from; i >= to; i-- ) rx[i] = xk;
    from = to - 1;
  }
  free( idx );
} // monotone

// [[Rcpp::export]]
Rcpp::NumericVector monotone(Rcpp::NumericVector x, Rcpp::NumericVector w) {

  double* x_new = x.begin();
  double* w_new = w.begin();
  int n = x.length();

  monotoneC(&n, x_new, w_new);

  return x;


}
