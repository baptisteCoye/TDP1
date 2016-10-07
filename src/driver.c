#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "util.h"
#include "cblas.h"
#include "perf.h"

int main(int argc, char ** argv){
  
  const enum CBLAS_ORDER Order = CblasColMajor;

  int m = 10, n = 7;
  double * A;
  int lda = allouer_matrice(&A, m, m);
  int i, j;

  //Initialisation de la Matrice

  for (i = 0; i < m; ++i){
    for (j = 0; j < m; ++j){
      A[i+j*lda] = (double) j;
    }
  }

  //Affichage de la Matrice
  printf("#%f#\n", A[0]);
  affiche(m, m, A, lda, stdout);

  double * X, * Y;
  m = 50;

  perf_t begin, end;

  while(m < 1000000){

    int incx = allouer_vecteur(&X, m);
    int incy = allouer_vecteur(&Y, m);
    //Test du Produit Scalaire
    double dot = cblas_ddot(m, X, incx, Y, incy);

    for (i = 0; i < m; ++i){
      X[i] = i;
      Y[i] = i;
    }

    perf(&begin);
    cblas_ddot(m, X, incx, Y, incy);
    perf(&end);

    printf("perf avec m = %d :\n", m);
    perf_diff(&begin, &end);
    printf("   Temps d'execution : "); perf_printmicro(&end);
    printf("   MFlops/s : %f\n\n", perf_mflops(&end, (long) m));

    desallouer_vecteur(X);
    desallouer_vecteur(Y);

    m += m/4;
  }
  
  m = 10;
  double * B;
  int ldb = allouer_matrice(&B, m, m);

  for (i = 0; i < m*m; ++i){
    B[i] = 1;
  }

  affiche(m, m, B, ldb, stdout);
  double * C;
  int ldc = allouer_matrice(&C, m, m);

  cblas_dgemm_scalaire(m, m, A, lda, B, ldb, C, ldc);

  affiche(m, m, C, ldc, stdout);

  desallouer_matrice(A);
  desallouer_matrice(B);
  desallouer_matrice(C);
  
}
