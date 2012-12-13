#define L1_BLK_SIZE 2
#define MULTIPLE 2
#define L2_BLK_SIZE (L1_BLK_SIZE * MULTIPLE)
#define verbose 0
#define l2squared = L2_BLK_SIZE * L2_BLK_SIZE

/* HPC 2012 Project : Jacqueline Bush, Paul Torres */

/* Code Preforms Blocked Matrix Matrix Multiplication, Matrices do not need to be square */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "MatrixVector.h"
#include "Utilities.h"
#include "test.h"

/* Blocked Matrix Matrix Multiply Subfunctions */
/* We couldn't add these two functions to the header because they are static */

static void dgemm_vector_lowest(const double * restrict A, const double * restrict B, 
				double * restrict C, const int h_padding, 
				const int w_padding);
static void dgemm_vector_middle(const double*restrict A, const double*restrict B, 
				double*restrict C, const int iblock, 
				const int jblock, const int wA, const int hA);
void AltBlockMatrixByCol(const double * const inA, double * outA, 
			 const int blocksize, const int i, const int j, 
			 const int hA);
int AltMatrixVector_middle(const double * const A, const int hA, const int wA,
			   const double * const B, double *C);
int MatrixVectorMultiply(const double * const A, const int hA, const int wA, 
			 const double * const B, double *C){
  int hB = wA;
  /* Matrix Deminsions of Output matrix */
  int wC = 1;
  int hC = hA;
  
  CleanMatrix(C, hC, wC);
  
  /* Calculate the number of Blocks in the width of A */
  int wn_bloc = (wA + L2_BLK_SIZE - 1) / L2_BLK_SIZE;
  // Number of Blocks in the height (round up)
  int hn_bloc = (hC+L2_BLK_SIZE -1)/L2_BLK_SIZE; 
  // Number of Blocks in the widith of matrix A  (round up)
  
  int b = L2_BLK_SIZE; // Block Size
  
  /* Make static matrices */
  static __attribute__ ((aligned(16))) double a_block[L2_BLK_SIZE*L2_BLK_SIZE];
  static __attribute__ ((aligned(16))) double b_block[L2_BLK_SIZE];
  static __attribute__ ((aligned(16))) double c_block[L2_BLK_SIZE];
  
  /* Clean matrices */
  CleanMatrix(c_block, b, 1);
  CleanMatrix(a_block, b, b);
  CleanMatrix(b_block, b, 1);
  
  for(int i = 0; i < hn_bloc; i++){
    for(int j = 0; j < wn_bloc; j++){
      BlockVectorForMatrixVector(C, c_block, hC, b, i);
      BlockMatrix(A, a_block, hA, wA, b, i, j);
      BlockVectorForMatrixVector(B, b_block, hB, b, j);
      if(verbose) printf("\n=top=i:%d j:%d a_block: \n",i,j);
      prettyPrint(a_block, b, b);
      if(verbose) printf("\n=top=i:%d j:%d b_block: \n", i,j);
      prettyPrint(b_block, b, 1);
      dgemm_vector_middle(a_block, b_block, c_block, i, j, wA, hA);
      CleanMatrix(a_block, b, b);
      CleanMatrix(b_block, b, 1);
      UnBlockVectorForMatrixVector(C, c_block, hC, b, i);
      if(verbose) printf("\n=top=C:\n");
      prettyPrint(C, wA, 1);
      if(verbose) printf("\n=top= c_block:\n");
      prettyPrint(c_block, b, 1);
      CleanMatrix(c_block, b, b);
    }
  }
  
  /* Uncomment to Test this function */
  //testMatrixMultiply(A, B, C, hA, wA, hB, 1);
  return 0;
}


/*=============================================================================*/
 

/*----------------------------dgemm_middle Code-------------------------------*/
static void dgemm_vector_middle(const double*restrict A, const double*restrict B, 
				double*restrict C, const int iblock, 
				const int jblock, const int wA, const int hA){
  /*----------------------------------------------------------------------------- 
    PURPOSE: Takes fixed sized matrices A and B, size L2_BLK_SIZE, Outputs result matrix C.  
    ARGUMENTS:
    -----------------------------------------------------------------------------*/
  int n = L2_BLK_SIZE; // Size of Matrix 
  int n_bloc = L2_BLK_SIZE/L1_BLK_SIZE; // Number of Blocks
  int b = L1_BLK_SIZE; // Block Size
  int h_padding = 0;
  int hn_blockl2 = (hA + L2_BLK_SIZE - 1) / L2_BLK_SIZE;
  //if we're in last L2 block along i of A cut short i along answer vector
  if(iblock == hn_blockl2){
    //figure out how much padding will be in last L1 block and pass it on to lowest
    h_padding = iblock * L2_BLK_SIZE - hA;
  }
  //if we're in last L2 block along j of A cut short i along input vector
  int w_padding = 0;
  int wn_blockl2 = (wA + L2_BLK_SIZE - 1) / L2_BLK_SIZE;
  if(jblock == wn_blockl2){
    w_padding = jblock * L2_BLK_SIZE - wA;
  }
  
  /* Make static matrices */
  static __attribute__ ((aligned(16))) double a_block[L1_BLK_SIZE*L1_BLK_SIZE];
  static __attribute__ ((aligned(16))) double b_block[L1_BLK_SIZE];
  static __attribute__ ((aligned(16))) double c_block[L1_BLK_SIZE]; 
 
  CleanMatrix(c_block, b, 1);
  for(int i = 0; i < n_bloc; i++){
    for(int j = 0; j < n_bloc; j++){
      BlockVectorForMatrixVector(C, c_block, n, b, i);
      BlockMatrix(A, a_block, n, n, b, i, j);
      BlockVectorForMatrixVector(B, b_block, n, b, j);
      if(verbose) printf("\n=middle=i:%d j:%d a_block: \n",i,j);
      prettyPrint(a_block, b, b);
      if(verbose) printf("\n=middle=i:%d j:%d b_block: \n", i,j);
      prettyPrint(b_block, b, 1);
      if(i == (n_bloc-1) && j == (n_bloc-1)){
	dgemm_vector_lowest(a_block, b_block, c_block, h_padding, w_padding);
      } else if(i == (n_bloc-1)){
	dgemm_vector_lowest(a_block, b_block, c_block, h_padding, 0);
      } else if(j == (n_bloc-1)){
	dgemm_vector_lowest(a_block, b_block, c_block, 0, w_padding);
      } else if(i < n_bloc && j < n_bloc){
	dgemm_vector_lowest(a_block, b_block, c_block, 0, 0);
      }
      CleanMatrix(a_block, b, b);
      CleanMatrix(b_block, b, 1);
      if(verbose) printf("\n=middle=i:%d j:%d C: \n", i, j);
      prettyPrint(C, n, n);
      UnBlockVectorForMatrixVector(C, c_block, n, b, i);
      if(verbose) printf("\n=middle=i:%d j:%d C: \n", i, j);
      prettyPrint(C, n, n);
      if(verbose) printf("\n=middle=i:%d j:%d c_block: \n", i, j);
      prettyPrint(c_block, b, b);
      CleanMatrix(c_block, b, 1);
    }
  }
}


/*----------------------------dgemm_lowest Code-------------------------------*/
/*----------------------------------------------------------------------------- 
  PURPOSE: Takes fixed sized matrix A and vector B, size L1_BLK_SIZE, Outputs 
  result vector C...avoids operations on zeroes
  -----------------------------------------------------------------------------*/
static void dgemm_vector_lowest(const double * restrict A, const double * restrict B, 
				double * restrict C, const int h_padding, 
				const int w_padding){
  int n = L1_BLK_SIZE - h_padding;
  int n2 = L1_BLK_SIZE - w_padding;
  if(verbose){
    printf("\n=Vector lowest=");
    printf("h_padding: %d, w_padding: %d", h_padding, w_padding);
    printf("\nA:\n");
    prettyPrint(A, n, n);
    printf("\nB:\n");
    prettyPrint(B, n, 1);
    printf("\nC:\n");
    prettyPrint(C, n, 1);
  }
  
  for(int j = 0; j < n2; j++){
    for(int i = 0; i < n ; i++){
      if(verbose) printf("\nj:%d i:%d  C[%d]: %f  C[i]= %f  A[i + j*n]: %f x B[j]: %f", j, i, i, C[i], A[i + j*n]*B[j], A[i + j*n], B[j]);
      C[i] += A[i + j*n] * B[j];
    }
  }
}

/*
  This accesses the elements of M sequentially unlike blocking version. Might be faster.
  A  -> Matrix that is m x n
  B  -> Vector m x 1
  hA -> m of A
  wA -> n of A
  C  -> the output vector
*/
/*
int AltMatrixVectorMultiply(const double * const A, const int hA, const int wA, 
			 const double * const B, double *C){
  static __attribute__ ((aligned(16))) double a_block[3 * L2_BLK_SIZE * L2_BLK_SIZE / 2];

  int bigblocksize = 3 * L2_BLK_SIZE  * L2_BLK_SIZE / 2;
  //compute # of L2 blocks for width and height of A
  int wn_blocA = (wA + bigblocksize - 1) / bigblocksize;
  //by definition
  int hn_blocB = wn_blocA;
  // Number of Blocks in the height
  int hn_blocA = (hA + bigblocksize - 1) / bigblocksize;

  for(int j = 0; j < wn_blocA; j++){
    for(int i = 0; i < hn_blocA; i++){
      //C[i] = A[i + j * hA] * B[i];
      AltBlockMatrixByCol(A, a_block, bigblocksize, i, j, hA);
      AltMatrixVectorMiddle(a_block, hA, wA, B[j], C);
      
    }
  }

}


void AltBlockMatrixByCol(const double * const inA, double * outA, 
			 const int blocksize, const int i, const int j, 
			 const int hA){
  int num_bloc = 0;
  int padding = 0;
  int loopguard = blocksize;
  num_bloc = (hA + blocksize - 1) / blocksize;
  if(i == num_bloc - 1){
    padding = num_bloc * blocksize - hA;
    loopguard = blocksize - padding;
  }
  for(int i2 = 0; i2 < loopguard; i2++){
    outA[i] = inA[i2 + i * blocksize + j * hA];
  }
}

int AltMatrixVector_middle(const double * const A, const int hA, const int wA,
			   const double B, double *C){
  static __attribute__ ((aligned(16))) double a_block1[3 * L1_BLK_SIZE * L1_BLK_SIZE / 2];
  int bigblocksize = 3 * L1_BLK_SIZE  * L1_BLK_SIZE / 2;
  //compute # of L1 blocks for width and height of A
  //int wn_blocA = (wA + bigblocksize - 1) / bigblocksize;
  //by definition
  int hn_blocB = wn_blocA;
  // Number of Blocks in the height
  int hn_blocA = (hA + bigblocksize - 1) / bigblocksize;

  for(int j = 0; j < wn_blocA; j++){
    for(int i = 0; i < hn_blocA; i++){
      //C[i] = A[i + j * hA] * B[i];
      AltBlockMatrixByCol(A, a_block1, bigblocksize, i, j, hA);
      
    }
  }

}
*/
