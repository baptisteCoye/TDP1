#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "util.h"
#include "cblas.h"
#include "perf.h"

int main(int argc, char ** argv){
  
  int m = 10, n = 7;
  double * A;
  int lda = allouer_matrice(&A, m, n);
  int i, j;

  //Initialisation de la Matrice

  for (i = 0; i < m; ++i){
    for (j = 0; j < n; ++j){
      A[i+j*lda] = (double) i;
    }
  }

  //Affichage de la Matrice
  affiche(m, n, A, lda, stdout);

  double * X, * Y;
  m = 50;

  perf_t begin, end;

  while(m < 1000000){
  //Initialisation des vecteurs
  for (i = 0; i < m; ++i){
    X[i] = i;
    Y[i] = i;
  }

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
  


  desallouer_matrice(A);
}
