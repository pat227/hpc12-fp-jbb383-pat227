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



