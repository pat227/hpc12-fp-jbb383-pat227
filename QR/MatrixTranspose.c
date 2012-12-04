#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "QR.h"

#define L1_BLK_SIZE 2
#define MULTIPLE 2
#define L2_BLK_SIZE (L1_BLK_SIZE * MULTIPLE)
void copyTransposedL2Block(const double * const B, const int submatrix_row , const int submatrix_col, double * A, const int w, const int h);
void copyL2Block(const double * const A, const int submatrix_row , const int submatrix_col, double * B, const int w, const int h);
void copyL1Block(const double * const A, const int i, const int j, double * B);
void transposeL2(const double * const A, double * B);
void copyTransposedL1Block(const double * const B, const int submatrix_row, const int submatrix_col, double * A);
void transposeL1Size(double * const A, double * B);
void prettyPrint(const double * const A, const int size){
  printf("\n");
  for(int i = 0; i < size ; i++){
    for(int j = 0; j < size ; j++){
      printf("%5.5f ", A[i+j*size]);
    }
    printf("\n");
  }
  printf("\n");
}
/*
  Transposes a matrix using L1 and L2 sized sub-blocks to hopefully gain a speed
  up, similar to the openCL matrix transposition although that was on a GPU...
*/
//B must have the dimensions of A transposed--I do not check for it; function transposes 
//the matrix by blocks; each element of each submatrix block is transposed and then each block
//is copied back to its transposed location within the larger matrix
void MatrixTranspose(const double * const A, const int w, const int h, double * B){
  int num_blocks_w = (w + L2_BLK_SIZE - 1) / L2_BLK_SIZE;
  int num_blocks_h = (h + L2_BLK_SIZE - 1) / L2_BLK_SIZE;
  printf("\n# of sub-blocks: %d by %d ", num_blocks_w, num_blocks_h);
  /* Make static matrices */
  static __attribute__ ((aligned(16))) double b_blockL2[L2_BLK_SIZE*L2_BLK_SIZE];
  static __attribute__ ((aligned(16))) double c_blockL2[L2_BLK_SIZE*L2_BLK_SIZE];

  for(int j = 0; j < num_blocks_w; j++){
    for(int i = 0; i < num_blocks_h; i++){
      copyL2Block(A, i , j, b_blockL2, w, h);
      printf("\nb_blockL2:\n");
      prettyPrint(b_blockL2, L2_BLK_SIZE);
      printf("\nInvoking tranposeL2");
      transposeL2(b_blockL2, c_blockL2);
      printf("\nc_blockL2:\n");
      prettyPrint(c_blockL2, L2_BLK_SIZE);
      copyTransposedL2Block(c_blockL2, i, j, B, w, h);
    }
  }
}

//copy a transposed L2-sized block to proper location in larger answer matrix--submatrices do not go
//back to where they came from but are themselves transposed as well
void copyTransposedL2Block(const double * const B, const int submatrix_row , const int submatrix_col, double * A, const int w, const int h){
  for(int j = 0; j < L2_BLK_SIZE; j++){
    for(int i = 0; i < L2_BLK_SIZE; i++){
      A[(submatrix_col*L2_BLK_SIZE) + i + (submatrix_row+j)*h] = B[i+j*L2_BLK_SIZE];
    }
  }
}
//copy a NOT-YET-transposed L2 submatrix block of A to L2-block-sized B
void copyL2Block(const double * const A, const int submatrix_row , const int submatrix_col, double * B, const int w, const int h){
  for(int j = 0; j < L2_BLK_SIZE ; j++){
    for(int i = 0; i < L2_BLK_SIZE ; i++){
      B[i+j*L2_BLK_SIZE] = A[(submatrix_row*L2_BLK_SIZE) + i + (submatrix_col+j)*h];
    }
  }
}
//copy a NOT-YET-transposed L1 submatrix block of L2-sized A to L1-block-sized B
void copyL1Block(const double * const A, const int i, const int j, double * B){
  printf("=CopyL1Block=");
  for(int j2 = 0; j2 < L1_BLK_SIZE ; j2++){
    for(int i2 = 0; i2 < L1_BLK_SIZE ; i2++){
      B[i2 + j2 * L1_BLK_SIZE] = A[i*L1_BLK_SIZE + (j*L1_BLK_SIZE)*L2_BLK_SIZE + i2 + j2*L2_BLK_SIZE];
    }
  }
}
//given an L2 sized submatrix A, transpose it, one L1 submatrix at a time, into B
void transposeL2(const double * const A, double * B){
  static __attribute__ ((aligned(16))) double b_blockL1[L1_BLK_SIZE*L1_BLK_SIZE];
  static __attribute__ ((aligned(16))) double c_blockL1[L1_BLK_SIZE*L1_BLK_SIZE]; 
  printf("\n=transposeL2=");
  for(int j = 0; j < MULTIPLE ; j++){
    for(int i = 0; i < MULTIPLE ; i++){
      copyL1Block(A, i ,j , b_blockL1);
      printf("\nb_blockL1:\n");
      prettyPrint(b_blockL1, L1_BLK_SIZE);
      transposeL1Size(b_blockL1, c_blockL1);
      printf("\nc_blockL1:\n");
      prettyPrint(c_blockL1, L1_BLK_SIZE);
      copyTransposedL1Block(c_blockL1, i, j, B);
    }
  }
}
//copy a tranposed L1-block sized block B into A which is L2 sized, while transposing submatrix itself
//amongst the other submatrices
void copyTransposedL1Block(const double * const B, const int submatrix_row, const int submatrix_col, double * A){
  for(int j = 0; j < L1_BLK_SIZE; j++){
    for(int i = 0; i < L1_BLK_SIZE; i++){
      A[(submatrix_row * L1_BLK_SIZE) + (submatrix_col*L2_BLK_SIZE) + i + j] = B[i+j*L1_BLK_SIZE];
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
