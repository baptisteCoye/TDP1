#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "util_float.h"
#include "cblas.h"
#include "perf.h"
#include "dgemm.h"

int main(int argc, char **argv){
  const enum CBLAS_ORDER Order = CblasColMajor;

  int m = 10, n = 7;
  float *A;
  int lda = allouer_matrice(&A, m, m);
  int i, j, k;

  //Initialisation de la Matrice

  for (i = 0; i < m; ++i){
    for (j = 0; j < m; ++j){
      A[i+j*lda] = (float) j;
    }
  }

  //Affichage de la Matrice
  printf("#%f#\n", A[0]);
  affiche(m, m, A, lda, stdout);

  float *X, *Y;
  m = 50;
  //TO DO: Test l'implémentation de Saxpy
  return 1;
}
