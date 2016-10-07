#include "cblas.h"

double cblas_ddot(const int m, const double * DX, const int INCX, 
		  const double * DY, const int INCY){
  double res = 0;
  int i, j;

  for (i = 0; i < m; ++i){
    res += DX[i] * DY[i];
  }

  return res;
}
