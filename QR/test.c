
/* HPC 2012 Final Project: Jacqueline Bush, Paul Torres */

/* This file contains code used to check rotines used in the project, such as Matrix Matrix Multiply and Matrix Transpose. We use the well know sequential routines to check the optimal routines we wrote. */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "MatrixMatrixMultiply.h"
#include "MatrixTranspose.h"
#include "test.h"



/*=================== Code Tests Matrix Matrix Multiply =======================*/


 /*----------------------------dgemm_simple Code-------------------------------*/
void dgemm_simple(const double *A, const int wA, const int hA, const double *B, const int wB, const int hB, double *C) {
/*----------------------------------------------------------------------------- 
PURPOSE: Computes simple matrix multiplication with A and B in Column-major order. 
ARGUEMENTS:
	wA: Width of A, number of columns in A
	hA: Height of A, number of rows in A
	wB: Width of B
	hB: Height of B 
-----------------------------------------------------------------------------*/
 
int hC = hA;
int wC = wB;

for(int i=0; i<hC; i++){
	for(int j=0;j<wC; j++){
		for(int k=0; k<wA;k++){
			C[i + j*hC] += A[i +k*hA]*B[k +j*hB];
		}
	}
  } 
}

void testMatrixMultiply( const double *A, const double *B, double *C, int hA, int wA, int hB, int wB ){
/*----------------------------------------------------------------------------- 
PURPOSE: Tests Blocked matrix multiplication against the simple variant (which we know works)
ARGUEMENTS:
	hA: Height of A
	wA: Width of A
	hB: Height of B
	wB: Width of B
-----------------------------------------------------------------------------*/

/* Dimensions of C */
int hC = hA;
int wC = wB;

/* Test Matrix */
double *testC = malloc( wC*hC*sizeof(double) );
CleanMatrix(testC, hC, wC);
dgemm_simple( A, wA, hA, B, wB, hB, testC);

		
/* Error Check */
for(int i=0; i<hC; i++){
	for(int j=0;j<wC; j++){
		double error = abs( C[i+ j*hC] - testC[i + j*hC]);
		//printf(" error = %f, i = %d, j = %d\n", error, i, j);
		//printf(" C = %f, testC = %f \n", C[i+j*hC], testC[i + j*hC]);
		double errorbound = 1e-5;
		if( error > errorbound ){
		fprintf(stderr,"Blocked Matrix Multiplication is not working! \n") ;
   	 	abort();
		}
	}
  }

}

/*============================================================================*/


