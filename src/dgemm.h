#ifndef DGEMM_H
#define DGEMM_H

#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "cblas.h"

/*void cblas_dgemm_scalaire(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
			  const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
			  const int K, const double alpha, const double *A,
			  const int lda, const double *B, const int ldb,
			  const double beta, double *C, const int ldc);
*/

void cblas_dgemm_scalaire(const int M,
			  const int N,
			  const double *A,
			  const int lda,
			  const double *B,
			  const int ldb,
			  double *C,
			  const int ldc);

void matmultKIJ( const int I, const int J, const int K,
		 const double *A, const int lda, 
		 const double *B, const int ldb, 
		 double *C, const int ldc);

void matmultIJK( const int I, const int J, const int K,
		 const double *A, const int lda, 
		 const double *B, const int ldb, 
		 double *C, const int ldc);

void matmultJIK( const int I, const int J, const int K,
		 const double *A, const int lda, 
		 const double *B, const int ldb, 
		 double *C, const int ldc);

#endif
