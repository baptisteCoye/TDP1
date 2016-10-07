#include "ddot.h"

double ddot(int m, double * DX, int INCX, double * DY, int INCY){
  double res = 0;
  int i, j;

  for (i = 0; i < m; ++i){
    res += DX[i] * DY[i];
  }

  return res;
}
