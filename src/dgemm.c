#include "dgemm.h"

void cblas_dgemm_scalaire( const int M, const int N,
			   const double *A, const int lda, 
			   const double *B, const int ldb,
			   double *C, const int ldc){



/* void cblas_dgemm_scalaire(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
			  const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
			  const int K, const double alpha, const double *A,
			  const int lda, const double *B, const int ldb,
			  const double beta, double *C, const int ldc){*/
  int i, j;

  printf("Call to cblas_dgemm_scalaire, variables are :\n");
  printf("M = %d\n", M);
  printf("N = %d\n", N);
  printf("lda = %d\n", lda);
  printf("ldb = %d\n", ldb);
  printf("ldc = %d\n", ldc);

  for(i = 0; i < M; ++i){
    for(j = 0; j < M; ++j){
      C[i*ldc+j] = cblas_ddot(M, &(A[j*lda]), lda, &(B[i*ldb]), ldb);
    }
  }
}


