#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "util.h"
#include "cblas.h"
#include "perf.h"
#include "dgemm.h"

int main(int argc, char ** argv){
  
  const enum CBLAS_ORDER Order = CblasColMajor;

  int m = 10, n = 7;
  double * A;
  int lda = allouer_matrice(&A, m, m);
  int i, j, k;

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
    cblas_daxpy(m, 2, X, incx, Y, incy);
    perf(&end);
    printf("perf daxpy avec m = %d :\n",m);
    perf_diff(&begin, &end);
    printf("   Temps d'execution : "); perf_printmicro(&end);
    printf("   MFlops/s : %f\n\n", perf_mflops(&end, (long) m));
    
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
  
  ////////////////////////////////////////////////////////////
  ///         Test des multiplications de matrices         ///
  ////////////////////////////////////////////////////////////

  printf("###########################################\n                 matmult\n###########################################\n\n");

  int M = 1000;

  perf_t kij, ijk, jik;
  perf_t totalkij, totalijk, totaljik;

#ifdef _OPENMP
  printf("Open mp ok\n");
#else 
  printf("Open mp ko\n");
#endif
  

  while (M == 1000){
    printf("##############\nperf avec M = %d :\n\n", M);
    lda = allouer_matrice(&A, M, M);
    ldb = allouer_matrice(&B, M, M);
    ldc = allouer_matrice(&C, M, M);

    for (i = 0; i < M; ++i){
      for (k = 0; k < M; ++k){
	A[i+k*lda] = (double) k;
      }
    }

    for (k = 0; k < M; ++k)
      for (j = 0; j < M; ++j)
	B[k+j*ldb] = 1;
 
    /* printf("A = \n"); */
    /* affiche(M,M,A,lda,stdout); */
    /* printf("B = \n"); */
    /* affiche(M,M,B,ldb,stdout); */


    perf(&begin);
    matmultKIJ(M,M,M, 
	       A, lda,
	       B, ldb,
	       C, ldc);
    perf(&kij);

    /* printf("Ckij = \n"); */
    /* affiche(M,M,C,ldc,stdout); */

    printf("KIJ : \n\n");
    perf_diff(&begin, &kij);
    printf("   Temps d'execution : "); perf_printmicro(&kij);
    printf("   MFlops/s : %f\n\n", perf_mflops(&kij, (long) M*M*M));

    perf(&begin);
    matmultIJK(M,M,M, 
	       A, lda,
	       B, ldb,
	       C, ldc);
    perf(&ijk);
    /* printf("Cijk = \n"); */
    /* affiche(M,M,C,ldc,stdout); */
    
    printf("IJK : \n\n");
    printf("perf avec M = %d :\n", M);
    perf_diff(&begin, &ijk);
    printf("   Temps d'execution : "); perf_printmicro(&ijk);
    printf("   MFlops/s : %f\n\n", perf_mflops(&ijk, (long) M*M*M));

    perf(&begin);
    matmultJIK(M,M,M, 
	       A, lda,
	       B, ldb,
	       C, ldc);
    perf(&jik);

    /* printf("Cjik = \n"); */
    /* affiche(M,M,C,ldc,stdout); */
    
    printf("JIK : \n\n");
    printf("perf avec M = %d :\n", M);
    perf_diff(&begin, &jik);
    printf("   Temps d'execution : "); perf_printmicro(&jik);
    printf("   MFlops/s : %f\n\n", perf_mflops(&jik, (long) M*M*M));

    desallouer_matrice(A);
    desallouer_matrice(B);
    desallouer_matrice(C);

    M += 100;
  }
  
}
