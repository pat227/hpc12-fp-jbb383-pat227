
/* HPC 2012 Project  : Jacqueline Bush, Paul Torres */

/* BlockedQR FILE: 
Performs Blocked QR factorization using static Matrices
*/

/* Headers: */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "WY.h"
#include "test.h"
#include "Utilities.h"
#include "BlockedQR.h"
#include "BlockedQR2.h"

#define L1_BLK_SIZE 16
#define L2_BLK_SIZE (L1_BLK_SIZE * 32)

static void dgemm_lowest( const double*restrict, const double*restrict, double*restrict);
static void dgemm_middle( const double*restrict, const double*restrict, double*restrict);
static void staticWY( const double*restrict A, double*restrict Q, double*restrict Qt, double*restrict R);
static void BlockedQRMiddle(const double*restrict A, double*restrict Q, double*restrict Qt);


void BlockedQR2( double *A, int h, int w, double *Q){

/* Set Q = m by m identity matrix */
 	for(int i = 0; i<h; i++){
		for(int j=0; j<h; j++){
		Q[i + j*h] = 0;
		if(i == j)
			Q[i +j*h] = 1;
		}
	}	





/* Make static matrices */
static __attribute__ ((aligned(16))) double a_block[L1_BLK_SIZE*L1_BLK_SIZE];
static __attribute__ ((aligned(16))) double b_block[L1_BLK_SIZE*L1_BLK_SIZE];
static __attribute__ ((aligned(16))) double c_block[L1_BLK_SIZE*L1_BLK_SIZE];



}

static void BlockedQRMiddle(const double*restrict A, double*restrict Q, double*restrict Qt)
{





}



/*=========================================================================*/

static void staticWY( const double*restrict A, double*restrict Q, double*restrict Qt, double*restrict R){





}


/*=========================================================================*/


/* Matrix Multiplication Code */ 
/*----------------------------dgemm_middle Code-------------------------------*/
static void dgemm_middle(const double*restrict A, const double*restrict B, double*restrict C){
/*----------------------------------------------------------------------------- 
PURPOSE: Takes fixed sized matrices A and B, size L2_BLK_SIZE, Outputs result matrix C.  
ARGUEMENTS:
-----------------------------------------------------------------------------*/

int n = L2_BLK_SIZE; // Size of Matrix 
int n_bloc = L2_BLK_SIZE/L1_BLK_SIZE; // Number of Blocks
int b = L1_BLK_SIZE; // Block Size

/* Make static matrices */
static __attribute__ ((aligned(16))) double a_block[L1_BLK_SIZE*L1_BLK_SIZE];
static __attribute__ ((aligned(16))) double b_block[L1_BLK_SIZE*L1_BLK_SIZE];
static __attribute__ ((aligned(16))) double c_block[L1_BLK_SIZE*L1_BLK_SIZE]; 

CleanMatrix(c_block, b, b);

for( int j = 0; j < n_bloc; j++){
	for( int i=0; i < n_bloc; i++){
		BlockMatrix(C, c_block, n, n, b, i, j);
		for( int k = 0; k < n_bloc ; k++){
		BlockMatrix(A, a_block, n, n, b, i, k);
		BlockMatrix(B, b_block, n, n, b, k, j);
		dgemm_lowest( a_block, b_block, c_block);
		CleanMatrix(a_block, b, b);
		CleanMatrix(b_block, b, b);	
		}
		UnBlockMatrix(C, c_block, n, n, b, i,j);		
	}
}	

}

/*----------------------------dgemm_lowest Code-------------------------------*/
static void dgemm_lowest(const double*restrict A, const double*restrict B, double*restrict C){
/*----------------------------------------------------------------------------- 
PURPOSE: Takes fixed sized matrices A and B, size L1_BLK_SIZE, Outputs result matrix C.  
ARGUEMENTS:
-----------------------------------------------------------------------------*/

int n = L1_BLK_SIZE , i, j, k;
 
  for(j=0; j<n; j++){
	for(k=0;k<n; k++){
		for(i=0; i<n;i++){
		C[i + j*n] += A[i + k*n]*B[k + j*n];
		}
	}
  }

}

/*============================================================================*/








