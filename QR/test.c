
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
#define EPSILON 0.00001

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

/*
  Function fails fast if any error occurs and doesn't waste time computing rest 
  of the matrix. Function does not store the result of the multiplication, but 
  only compares individually computed elements to elements in A. To support 
  non-square R, require more args to specify the dimensions of Q and R each.
  Perhaps structs with bounds would be better and they would know how to multiply
  against each other... 
  Q -> pointer to array of doubles in column-major order that represents Q
  R -> ditto that represents R
  A -> ditto that represents A
  Qm, Qn -> the m x n size of Q
  Rm, Rn -> the m x n size of R
*/
int IsQRequalToA(const double * const Q, const double * const R, 
		 const double * const A, const int Qm, const int Qn, 
		 const int Rm, const int Rn){
  //for each row i and col j in result
  int verbose = 0;
  double temp = 0.0;
  for(int i=0; i<Qm; i++){
    for(int  j=0; j<Rn; j++){
      for(int k=0; k<Qn; k++){
	temp += Q[i+k*Qm] * R[k+j*Rm];
	if(verbose) printf("\ni:%d j:%d k:%d +=  %-9.9f x %-9.9f",
			   i,j,k,Q[i+k*Qm],R[k+j*Rm]);
      }
      if(fabs(temp - A[i+j*Qm]) > EPSILON){
	if(verbose) printf("\nA != QR; element i: %d j: %d should be %-9.15f but computed to be %-9.15f\n",
	       i,j,A[i+j*Qm],temp);
	return 0;
      } else {
	temp = 0.0;
      }
    }
  }
  return 1;
}


/*============================================================================*/
