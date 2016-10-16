#include <assert.h>

#include "cblas.h"
#include "dgemm.h"

double cblas_ddot(const int m, const double * DX, const int INCX, 
		  const double * DY, const int INCY){
  double res = 0;
  int i, j;

  
  for (i = 0; i < m; ++i){
    res += DX[i] * DY[i];
  }

  return res;
}

void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc){

  assert((Order == CblasColMajor) && (TransB == CblasNoTrans));

  int i, j, k;

  

}

void cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY){
#pragma omp parallel for default(none) shared(X, Y)
  for (int i = 0; i < N; ++i){
    Y[i*incY] += alpha*X[i*incX];
  }
}

void cblas_saxpy(const int N, const float alpha, const float *X,
                 const int incX, float *Y, const int incY){
#pragma omp parallel for default(none) shared(X, Y)
  for (int i = 0; i < N; ++i){
    Y[i*incY] += alpha*X[i/incX];
  }
}

/*void cblas_dgemm_scalaire(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
			  const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
			  const int K, const double alpha, const double *A,
			  const int lda, const double *B, const int ldb,
			  const double beta, double *C, const int ldc){
*/
void cblas_dgemm_scalaire(const int M,
			  const int N,
			  const double *A,
			  const int lda,
			  const double *B,
			  const int ldb,
			  double *C,
			  const int ldc){
  int i, j;

  printf("Call to cblas_dgemm_scalaire, variables are :\n");
  printf("M = %d\n", M);
  printf("N = %d\n", N);
  printf("lda = %d\n", lda);
  printf("ldb = %d\n", ldb);
  printf("ldc = %d\n", ldc);

  for(i = 0; i < N; ++i){
    for(j = 0; j < M; ++j){
      C[i*ldc+j] = cblas_ddot(M, &(A[j*lda]), lda, &(B[i*ldb]), ldb);
    }
  }
}

void matmultKIJ( const int I, const int J, const int K,
		 const double *A, const int lda, 
		 const double *B, const int ldb, 
		 double *C, const int ldc){
  int i, j, k;

#pragma omp parallel for default(none) private(i,j) shared(A,B,C)
  for (j = 0; j < J; ++j){
    for (i = 0; i < I; ++i){
      C[i+j*ldc] = 0;
    }
  }


  for (k = 0; k < K; ++k){
#pragma omp parallel for default(none) private(i, j) shared(A,B,C, k)
    for (i = 0; i < I; ++i){
      for (j = 0; j < J; ++j){
	C[i+j*ldc] += A[i+k*lda]*B[k+ldb*j];
      }
    }
  }
}

void matmultIJK( const int I, const int J, const int K,
		 const double *A, const int lda, 
		 const double *B, const int ldb, 
		 double *C, const int ldc){

  int i, j, k;

#pragma omp parallel for default(none) private (i,j,k) shared(A,B,C)
  for (i = 0; i < I; ++i){
    for (j = 0; j < J; ++j){
      C[i+j*ldc] = 0;
      for (k = 0; k < K; ++k){
	C[i+j*ldc] += A[i+k*lda]*B[k+ldb*j];
      }
    }
  }  
}

void matmultJIK( const int I, const int J, const int K,
		 const double *A, const int lda, 
		 const double *B, const int ldb, 
		 double *C, const int ldc){
  int i, j, k;

#pragma omp parallel for default(none) private (i,j,k) shared(A,B,C)
  for (j = 0; j < J; ++j){
    for (i = 0; i < I; ++i){
      C[i+j*ldc] = 0;
      for (k = 0; k < K; ++k){
	C[i+j*ldc] += A[i+k*lda]*B[k+ldb*j];
      }
    }
  }  
}
