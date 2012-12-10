/* HPC 2012 Project : Jacqueline Bush, Paul Torres */

/* Code Preforms Blocked Matrix Matrix Multiplication, Matrices do not need to be square */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "MatrixMatrixMultiply.h"
#include "Utilities.h"
#include "test.h"

#define L1_BLK_SIZE 2
#define L2_BLK_SIZE (L1_BLK_SIZE * 2)

/* Blocked Matrix Matrix Multiply Subfunctions */
/* We couldn't add these two functions to the header because they are static */

static void dgemm_lowest( const double*restrict, const double*restrict, double*restrict);
static void dgemm_middle( const double*restrict, const double*restrict, double*restrict);

int MatrixMatrixMultiply( double *A, int hA, int wA, double *B, int hB, int wB, double *C)
{

/* Check that Matrix deminsions are valid */
if( wA !=  hB){
{
    fprintf(stderr, "Matrix Multiplication Invalid! \nThe number of columns of A must equal the number of rows of B!\n");
    abort(); 
  }
}


/* Matrix Deminsions of Output matrix */
int wC = wB;
int hC = hA;

CleanMatrix(C, hC, wC);

/* Calculate the number of Blocks in the height and width of C, as well as the width of A */
int wn_bloc = (wC+L2_BLK_SIZE -1)/L2_BLK_SIZE; // Number of Blocks in the widith (round up)
int hn_bloc = (hC+L2_BLK_SIZE -1)/L2_BLK_SIZE; // Number of Blocks in the height (round up)
int wA_bloc = (wA+L2_BLK_SIZE -1)/L2_BLK_SIZE; // Number of Blocks in the widith of matrix A  (round up)

int b = L2_BLK_SIZE; // Block Size

/* Make static matrices */
static __attribute__ ((aligned(16))) double a_block[L2_BLK_SIZE*L2_BLK_SIZE];
static __attribute__ ((aligned(16))) double b_block[L2_BLK_SIZE*L2_BLK_SIZE];
static __attribute__ ((aligned(16))) double c_block[L2_BLK_SIZE*L2_BLK_SIZE]; 

/* Clean matrices */
CleanMatrix(c_block, b, b);
CleanMatrix(a_block, b, b);
CleanMatrix(b_block, b, b);

for( int i = 0; i < hn_bloc; i++){
	for( int j=0; j < wn_bloc; j++){
		BlockMatrix(C, c_block, hC, wC,  b, i, j);
		for( int k = 0; k < wA_bloc ; k++){
		BlockMatrix(A, a_block, hA, wA, b, i, k);
		BlockMatrix(B, b_block, hB, wB, b, k, j);
		dgemm_middle(a_block, b_block, c_block);
		CleanMatrix(a_block, b, b);
		CleanMatrix(b_block, b, b);	
		}
		UnBlockMatrix(C, c_block, hC, wC, b, i,j);
		CleanMatrix(c_block, b, b);		
	}
}


/* Uncomment to Test this function */
testMatrixMultiply(A, B, C, hA, wA, hB, wB);


return 0;
}


/*=============================================================================*/
 

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




