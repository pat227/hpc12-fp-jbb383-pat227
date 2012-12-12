
/* HPC 2012 Final Project: Jacqueline Bush, Paul Torres */

/* This file contains code used to check rotines used in the project, such as Matrix Matrix Multiply and Matrix Transpose. We use the well know sequential routines to check the optimal routines we wrote. */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "MatrixMatrixMultiply.h"
#include "Utilities.h"
#include "MatrixTranspose.h"
#include "test.h"
#include "WY.h"


/*================= Code Test if Matrx is Upper Triangular ====================*/

void testUpperTriangular(double *R, int h, int w){

double errorbound = 1e-10;
int n = w;
if ( h < w)
	n= h;

  for(int i=0; i<n; i++){
    for(int j=0; j<i; j++){
      if(fabs(R[i+j*h]) > errorbound ){
	fprintf(stderr,"R is not uppertriangular! Fix error! \n") ;
   	 abort();
	}
   } }


}


/*================== Code Tests if Matrix is Orthogonal =======================*/

void testOrthogonal(double *Q, double *Qt, int h){
	
	double *testQ = malloc( h* h* sizeof(double));
	
	MatrixMatrixMultiply(Q, h,h, Qt, h, h, testQ);

	for(int i = 0; i<h; i++){
		for(int j =0 ; j< h; j++){
		double error = testQ[i + j*h] - 1;
		if (i !=j)
			error = testQ[i + j*h];
		double errorbound = 1e-5;
		if( error > errorbound ){
		fprintf(stderr,"Q is not orthogonal! Fix error! \n") ;
   	 	abort();
		}
	} }

}	 
  
/*=================== Code Tests Matrix Matrix Multiply =======================*/


 /*----------------------------dgemm_simple Code-------------------------------*/
void dgemm_simple(const double *A, const int hA, const int wA, const double *B, const int hB, const int wB, double *C) {
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

CleanMatrix(C , hC, wC);

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
dgemm_simple( A, hA, wA, B, hB, wB, testC);

		
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


/*====================== Code to test Matrix Transpose =======================*/

void simple_transpose(const double *A, int h, int w, double *B){
/*----------------------------------------------------------------------------- 
PURPOSE: Computes simple matrix transpose with A in Column-major order. 
ARGUEMENTS:
	w = width of A
	h = height of A
-----------------------------------------------------------------------------*/
  
  for(int j = 0; j < w; j++){
    for(int i = 0; i < h; i++){
      B[j + i*w] = A[i + j*h];
    }
  }
}


void testMatrixTranspose(const double *A, int h, int w, const double *ATranspose){
/*----------------------------------------------------------------------------- 
PURPOSE: Tests Blocked matrix transpose against the simple variant (which we know works)
ARGUEMENTS:
	h: Height of A
	w: Width of A
-----------------------------------------------------------------------------*/

/* Test Matrix */
double *testA = malloc( w*h*sizeof(double) );
CleanMatrix(testA, h, w);
simple_transpose( A, h,w, testA);

		
/* Error Check */
for(int i=0; i<w; i++){
	for(int j=0;j<h; j++){
		double error = abs( ATranspose[i+ j*w] - testA[i + j*w]);
		double errorbound = 1e-5;
		if( error > errorbound ){
		fprintf(stderr,"Blocked Matrix Tanspose is not working! \n") ;
   	 	abort();
		}
	}
  }

}


/*============================================================================*/







