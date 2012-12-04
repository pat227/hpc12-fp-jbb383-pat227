#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "QR.h"

#define L1_BLK_SIZE 2
#define MULTIPLE 2
#define L2_BLK_SIZE (L1_BLK_SIZE * MULTIPLE)
//B must have the dimensions of A transposed--I do not check for it; function transposes 
//the matrix by blocks; each element of each submatrix block is transposed and then each block
//is copied back to its transposed location within the larger matrix
void MatrixTranspose(const double * const A, const int w, const int h, double * B){
  int num_blocks_w = (w + L2_BLK_SIZE - 1) / L2_BLK_SIZE;
  int num_blocks_h = (h + L2_BLK_SIZE - 1) / L2_BLK_SIZE;
  
  /* Make static matrices */
  static __attribute__ ((aligned(16))) double b_blockL2[L2_BLK_SIZE*L2_BLK_SIZE];
  static __attribute__ ((aligned(16))) double c_blockL2[L2_BLK_SIZE*L2_BLK_SIZE];
  static __attribute__ ((aligned(16))) double b_blockL1[L1_BLK_SIZE*L1_BLK_SIZE];
  static __attribute__ ((aligned(16))) double c_blockL1[L1_BLK_SIZE*L1_BLK_SIZE]; 

//fill another supplied matrix with the tranpose of this one; for fixed size
//only so no need to worry about corner cases
  for(int j = 0; j < num_blocks_w; j++){
    for(int i = 0; i < num_blocks_h; i++){
      copyL2Block(A, i , j, b_blockL2, w, h);
      transposeL2(b_blockL2, c_blockL2);
    }
  }
}
//copy a transposed L2-sized block to proper location in larger answer matrix--submatrices do not go
//back to where they came from
void copyTransposedL2Block(const double * const A, const int submatrix_row , const int submatrix_col, double * B, const int w, const int h){
  for(int j = 0; j < L2_BLK_SIZE; j++){
    for(int i = 0; i < L2_BLK_SIZE; i++){
      A[(submatrix_col*L2_BLK_SIZE) + i + (submatrix_row+j)*h)] = B[i+j*L2_BLK_SIZE];
  }
}

//copy a NOT-YET-transposed L2 submatrix block of A to L2-block-sized B
void copyL2Block(const double * const A, const int submatrix_row , const int submatrix_col, double * B, const int w, const int h){
  for(int j = 0; j < L2_BLK_SIZE ; j++){
    for(int i = 0; i < L2_BLK_SIZE ; i++){
      B[i+j*L2_BLK_SIZE] = A[(submatrix_row*L2_BLK_SIZE) + i + (submatrix_col+j)*h)];
    }
  }
}

//copy a NOT-YET-transposed L1 submatrix block of L2-sized A to L1-block-sized B
void copyL1Block(const double * const A, double * B){
  for(int j = 0; j < MULTIPLE ; j++){
    for(int i = 0; i < MULTIPLE ; i++){
      for(int j2 = 0; j < L1_BLK_SIZE ; j++){
	for(int i2 = 0; i < L1_BLK_SIZE ; i++){
	  B[i2+j2*L1_BLK_SIZE] = A[];
	}
      }
    }
  }
}
//given an L2 sized submatrix, transpose it one L1 submatrix at a time
void transposeL2(const double * const A, double * B){
  for(int j = 0; j < MULTIPLE ; j++){
    for(int i = 0; i < MULTIPLE ; i++){  
      copyL1Block(A, i ,j , b_blockL1, const int w, const int h)
      transposeL1Size(b_blockL1, c_blockL1);
    }
  }
}

//fill another supplied matrix with the tranpose of this one; for fixed size
//only so no need to worry about corner cases
void transposeL1Size(double * const A, double * B){
  for(int j = 0; j < L1_BLK_SIZE; j++){
    for(int i = 0; i < L1_BLK_SIZE; i++){
      B[i+j*L1_BLK_SIZE] = A[j+i*L1_BLK_SIZE];
    }
  }
}
