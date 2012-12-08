/* HPC 2012 Final Project: Jacqueline Bush, Paul Torres */

/* Functions that didn't really have a place anywhere else were placed here. */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "Utilities.h"


void prettyPrint(const double * const A, const int m, const int n){
  for(int i = 0; i < m ; i++){
    for(int j = 0; j < n; j++){
      printf("%5.5f ", A[i+j*m]);
    }
    printf("\n");
  }
  printf("\n");
}

//h is height NOT how many; how many is start row - h, computed by function
void copyColVector(double * source, double * target, int start, int h, int col){
  for(int i = start ; i<h; i++){
    target[i + col*h] = source[i];
  }
}

//h is height NOT how many; how many is start row - h, computed by function
void copyMatrix(double * source, double * target, int h, int w){
  for(int j = 0 ; j<w; j++){
    for(int i = 0; i<h; i++){
      target[i + j*h] = source[i + j*h];
    }
  }
}
