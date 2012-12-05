/* HPC 2012 Project : Jacqueline Bush, Paul Torres */

/* Code Preforms Blocked Matrix Matrix Multiplication, Matrices do not need to be square */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "MatrixMatrixMultiply.h"

#define L1_BLK_SIZE 16
#define L2_BLK_SIZE (L1_BLK_SIZE * 32)


/* Blocked Matrix Matrix Multiply Subfunctions */
static void dgemm_lowest( const double*restrict, const double*restrict, double*restrict);
static void dgemm_middle( const double*restrict, const double*restrict, double*restrict);

int MatrixMatrixMultiply( double *A, int wA, int hA, double *B, int wB, int hB, double *C)
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
//test(A, B, C, hA, wA, hB, wB);

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

void test( const double *A, const double *B, double *C, int hA, int wA, int hB, int wB ){
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



/*---------------------Code to work with blocks of Matrix---------------------*/
 
void BlockMatrix(const double *inA, double *outA, int hA, int wA, int b, int i_bloc, int j_bloc){
/*----------------------------------------------------------------------------- 
PURPOSE: Takes matrix A, outputs the desired block, pads with zeros when necessary. 
ARGUEMENTS:
	hA: Height of matrix
	wA: Width of matrix  
	b: Blocksize
	i_bloc: i Block index 
	j_bloc: j Block index  
-----------------------------------------------------------------------------*/
int i, j;

if( wA % b != 0 || hA % b != 0  ){
	// Will only enter this section if we are in outerblock of code. 
	int wn_bloc = (wA+L2_BLK_SIZE -1)/L2_BLK_SIZE; // Number of Blocks 
	int hn_bloc = (hA+L2_BLK_SIZE -1)/L2_BLK_SIZE; // Number of Blocks
	int wpadding = wn_bloc*b - wA ; // Number of columns of zeros needed.
	int hpadding = hn_bloc*b - hA ; // Number of rows of zeros needed.
	for(i=0; i<b; i++){
		for(j=0; j<b; j++){
		
		if( (i_bloc == hn_bloc -1) && (i>= (b-hpadding))){
			outA[i+ j*b] = 0;
		} else{
			if( (j_bloc == wn_bloc -1) && (j >= (b-wpadding))){
			outA[i+ j*b] = 0;
			} else{
			outA[i + j*b] = inA[ i+i_bloc*b + (j+ j_bloc*b)*hA ] ;
			} 
		} }
    	} 
}else{
	for(i=0; i<b; i++){
		for(j=0; j<b; j++){
		outA[i + j*b] = inA[ i+i_bloc*b + (j+ j_bloc*b)*hA ] ;
		} 
	}
}

}


void UnBlockMatrix(double *outA, const double *inA, int hA, int wA, int b, int i_bloc, int j_bloc){
/*----------------------------------------------------------------------------- 
PURPOSE: Takes block of A, Replaces the values of A. Gets rid of pads if necessary.  
ARGUEMENTS:
	n: Size of matrix 
	b: Blocksize
	i_bloc: i Block index 
	j_bloc: j Block index  
-----------------------------------------------------------------------------*/

   int i, j;

        int wn_bloc = (wA+b -1)/b; // Number of Blocks 
	int hn_bloc = (hA+b -1)/b; // Number of Blocks
	int wpadding = wn_bloc*b - wA ; // Number of columns of zeros needed.
	int hpadding = hn_bloc*b - hA ; // Number of rows of zeros needed.

    if( ( wpadding == 0 && hpadding == 0) || (i_bloc != hn_bloc-1 && j_bloc != wn_bloc -1)){
	for(i=0; i<b; i++){
		for(j=0; j<b; j++){
		outA[ i+i_bloc*b + (j+j_bloc*b)*hA ] = inA[i + j*b];
		}
	} 
    }else{
	if( i_bloc == hn_bloc -1 && j_bloc != wn_bloc -1){
	for(i=0; i<(b-hpadding); i++){
		for(j=0; j<b; j++){
		outA[ i+i_bloc*b + (j+j_bloc*b)*hA ] = inA[i + j*b];
		}
	}
	}
	if( j_bloc == wn_bloc -1 && i_bloc != hn_bloc -1){
	for(i=0; i<b; i++){
		for(j=0; j<(b-wpadding); j++){
		outA[ i+i_bloc*b + (j+j_bloc*b)*hA ] = inA[i + j*b];
		}
	}
	}
	if( i_bloc == hn_bloc -1 && j_bloc == wn_bloc -1){
	for(i=0; i<(b-hpadding); i++){
		for(j=0; j<(b-wpadding); j++){
		outA[ i+i_bloc*b + (j+j_bloc*b)*hA ] = inA[i + j*b];
		}
	}
    	}
   }	
}

void CleanMatrix( double *A, int wA, int hA){
/*------------------------------------------------------------------------------
PURPOSE: Takes matrix A, Replaces it with zeros
ARGUEMENTS:
	wA: width of A
	hA: height of A
------------------------------------------------------------------------------*/

int i, j;

for(i=0; i<hA; i++){
	for(j=0; j<wA; j++){
	A[ i + j*hA] = 0 ;
	}
}

}



