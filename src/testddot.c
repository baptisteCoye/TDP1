#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "util.h"
#include "cblas.h"

int main(int argc, char ** argv){
  double * X, * Y;
  int i;
  int m = 50;

  while(m > 1000000){

    int incx = allouer_vecteur(&X, m);
    int incy = allouer_vecteur(&Y, m);

    for (i = 0; i < m; ++i){
      X[i] = i;
      Y[i] = i;
    }

    cblas_ddot(m, X, incx, Y, incy);

    desallouer_vecteur(X);
    desallouer_vecteur(Y);

    m += m/4;
  }

}
