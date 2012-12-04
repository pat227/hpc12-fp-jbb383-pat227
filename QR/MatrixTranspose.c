#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "QR.h"

#define L1_BLK_SIZE 2
#define L2_BLK_SIZE (L1_BLK_SIZE * 2)
//B must have the dimensions of A transposed--I do not check for it
void MatrixTranspose(double * A, int w, int h, double * B){
  int num_blocks_w = (w + L2_BLK_SIZE - 1) / L2_BLK_SIZE;
  int num_blocks_h = (h + L2_BLK_SIZE - 1) / L2_BLK_SIZE;
  
  

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
