#include <assert.h>
#include <math.h>
#include <omp.h>

#include "cblas.h"
#include "dgemm.h"

double cblas_ddot(const int m, const double * DX, const int INCX, 
		  const double * DY, const int INCY){
  double res = 0;
  int i, j;

  
  for (i = 0; i < m; ++i){
    res += DX[i*INCX] * DY[i*INCY];
  }

  return res;
}

void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc){

  assert(Order == CblasColMajor);

  int i, j, k;
  int Nt = 4;

  /*#ifdef _OPENMP
    Nt = omp_get_num_threads();
    #endif*/

  printf("Nt = %d\n", Nt);

  int sizeM[Nt], beginM[Nt];
  int sizeN[Nt], beginN[Nt];
  int sizeK[Nt], beginK[Nt];
  beginM[0] = 0;
  beginN[0] = 0;
  beginK[0] = 0;

  int tmpM = M;
  int tmpN = N;
  int tmpK = K;

  int tmpNt = Nt;

  for (i = 0; i < Nt; ++i){
    sizeM[i] = (tmpM + tmpNt - 1) / tmpNt;
    sizeN[i] = (tmpN + tmpNt - 1) / tmpNt;
    sizeK[i] = (tmpK + tmpNt - 1) / tmpNt;
    if (i > 0){
      beginM[i] = beginM[i-1] + sizeM[i-1];
      beginN[i] = beginN[i-1] + sizeN[i-1];
      beginK[i] = beginK[i-1] + sizeK[i-1];
    }
    tmpM -= sizeM[i];
    tmpN -= sizeN[i];
    tmpK -= sizeK[i];    
    tmpNt--;
  }

  printf("size M : ");
  for (i = 0; i < Nt; ++i){
    printf("%d ", sizeM[i]);
  }
  printf("\n");

  printf("begin M : ");
  for (i = 0; i < Nt; ++i){
    printf("%d ", beginM[i]);
  }
  printf("\n");

  for (i = 0; i < Nt; ++i){
    for (j = 0; j < Nt; ++j){
      for (k = 0; k < Nt; ++k){
	cblas_dgemm_scalaire(CblasColMajor, CblasNoTrans, CblasNoTrans,
			     sizeM[i], sizeN[j], sizeK[k],
			     1, &(A[beginK[k]*lda+beginM[i]]), lda, 
			     &(B[beginK[k]+beginN[j]*ldb]), ldb,
			     1, &(C[beginM[i]+beginN[j]*ldc]), ldc);
	printf("i = %d, sizeM[i] = %d\n", i, sizeM[i]);
      }
    }
  }
}


void cblas_dgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY){
  assert(Order == CblasColMajor);

  int i;

  if (TransA == CblasNoTrans){
#pragma omp parallel for default(none) private(i) shared(Y)
    for (i = 0; i < M; ++i){
      Y[i*incY] = beta*Y[i*incY] + alpha*cblas_ddot(N, &(A[i]), lda, X, incX);
    }
  } else if (TransA == CblasTrans) {
#prama omp parallel for default(none) private(i) shared(Y)
    for (i = 0; i < N; ++i)
      Y[i*incY] = beta*Y[i*incY] + alpha*cblas_ddot(M, &(A[i*lda], 1, X, incX);
  } else {
    fprintf(stderr, "cblas_dgemv erreur : CblasConjTrans non implémenté.\n");
  }
}

void cblas_dger(const enum CBLAS_ORDER order, const int M, const int N,
                const double alpha, const double *X, const int incX,
                const double *Y, const int incY, double *A, const int lda){
  assert(order == CblasColMajor);

  for (int i = 0; i < M; i++){
    for (int j = 0; j < N; j++){
      A[i+j*lda] += alpha*X[i*incX]*Y[j*incY];
    }
  }
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

void cblas_dgemm_scalaire(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
			  const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
			  const int K, const double alpha, const double *A,
			  const int lda, const double *B, const int ldb,
			  const double beta, double *C, const int ldc){
  /*
    void cblas_dgemm_scalaire(const int M,
    const int N,
    const double *A,
    const int lda,
    const double *B,
    const int ldb,
    double *C,
    const int ldc){*/

  assert(Order == CblasColMajor);
  
  int i, j;

  printf("Call to cblas_dgemm_scalaire, variables are :\n");
  printf("M = %d\n", M);
  printf("N = %d\n", N);
  printf("lda = %d\n", lda);
  printf("ldb = %d\n", ldb);
  printf("ldc = %d\n", ldc);

  if ((TransA == CblasNoTrans) && (TransB == CblasNoTrans)){
    for (j = 0; j < N; ++j){
      for (i = 0; i < M; ++i){
	C[i+j*ldc] = beta*C[i+j*ldc] + alpha*cblas_ddot(K, &(A[i]), lda, &(B[j*ldb]), 1); 
      }
    }
  } else if ((TransA == CblasTrans) && (TransB == CblasNoTrans)){
    for(j = 0; j < N; ++j){
      for(i = 0; i < M; ++i){
	C[i+j*ldc] = beta*C[i+j*ldc] + alpha*cblas_ddot(K, &(A[i*lda]), 1, &(B[i*ldb]), 1);
      }
    }    
  } else if ((TransA == CblasTrans) && (TransB == CblasTrans)){
    for (j = 0; j < N; ++j){
      for (i = 0; i < M; ++i){
	C[i+j*ldc] = beta*C[i+j*ldc] + alpha*cblas_ddot(K, &(A[i*lda]), 1, &(B[j]), ldb);
      }
    }
  } else if ((TransA == CblasNoTrans) && (TransB == CblasTrans)){
    for(j = 0; j < N; ++j){
      for(i = 0; i < M; ++i){
	C[i+j*ldc] = beta*C[i+j*ldc] + alpha*cblas_ddot(K, &(A[i]), lda, &(B[j]), ldb);
      }
    }
  } else {
    fprintf(stderr, "your choice of TransA and TransB is not implemented.\n");
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
