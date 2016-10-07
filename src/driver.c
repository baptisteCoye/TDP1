#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "util.h"
#include "ddot.h"

int main(int argc, char ** argv){
  
  int m = 10, n = 7;
  double * A;
  int lda = allouer_matrice(&A, m, n);
  int i, j;

  for (i = 0; i < m; ++i){
    for (j = 0; j < n; ++j){
      A[i+j*lda] = (double) i;
    }
  }

  affiche(m, n, A, lda, stdout);


  double * X, * Y;

  int incx = allouer_vecteur(&X, m);
  int incy = allouer_vecteur(&Y, m);

  for (i = 0; i < m; ++i){
    X[i] = i;
    Y[i] = i;
  }

  double dot = ddot(m, X, incx, Y, incy);

  printf ("    X = ");
  affiche(1, m, X, incx, stdout);
  printf ("    Y = ");
  affiche(1, m, Y, incy, stdout);

  printf("     X * Y = %10.0f\n", dot);

  desallouer_matrice(A);
}
