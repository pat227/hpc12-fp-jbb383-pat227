//#include<cstdlib.h> //rand, srand, malloc
#include<stdio.h> //fprintf, printf
#include<stdint.h> //exact types
#include<math.h>  //fabs
/*We use column major access and assume args are stored in column major order in accordance with all the other work we did this semester*/
//the rounding error tolerance
#define EPSILON 0.000001

//cannot be static when testing
/*inline static */ double identity(const uint32_t m, const uint32_t n){
  if(m==n){ 
    return 1.0;
  } else { 
    return 0.0;
  }
}
/*May have to use epsilon here since again getting two doubles to equal is sometimes impossible but they will be very nearly same value. Apply to Q in A=QR.*/
/*May have to use epsilon here since again getting two doubles to equal is sometimes impossible but they will be very nearly same value. Apply to Q of the A=QR.  Uses the smaller of the width or height. */
int IsQbyQtransposeIdentity(const double * const Q, const uint32_t m, const uint32_t n){
  /*  uint32_t n2 = 0;
  uint32_t stride = 0;
  if(m > n){
    n2 = n;
    stride = m;
  } else {
    n2 = m;
    stride = n;
    }*/
  for(uint32_t i=0; i<m; i++){
    for(uint32_t j=0; j<n; j++){
      for(uint32_t k=0; k<n; k++){
	if((Q[i+j*m] * Q[j+i*m]) != identity(i,j)){
	  printf("\nElement i: %d j: %d does not comport with I\n", i,j);
	  return 0;
	}
      }
    }
  }
  return 1;
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
int IsQRequalToA(const double * const Q, const double * const R, const double * const A, const uint32_t Qm, const uint32_t Qn, const uint32_t Rm, const uint32_t Rn){
  //for each row i and col j in result
  int verbose = 0;
  double temp = 0.0;
  for(uint32_t i=0; i<Qm; i++){
    for(uint32_t j=0; j<Rn; j++){
      for(uint32_t k=0; k<Qn; k++){
	temp += Q[i+k*Qm] * R[k+j*Rm];
	if(verbose) printf("\ni:%d j:%d k:%d +=  %-9.9f x %-9.9f",i,j,k,Q[i+k*Qm],R[k+j*Rm]);
      }
      if(fabs(temp - A[i+j*Qm]) > EPSILON){
	printf("\nA != QR; element i: %d j: %d should be %-9.15f but computed to be %-9.15f\n", i,j,A[i+j*Qm],temp);
	return 0;
      } else {
	temp = 0.0;
      }
    }
  }
  return 1;
}

/*Check constraint that a matrix is upper triangular. Apply to R only in A=QR.*/
int isUpperTriangular(const double * const M, const uint32_t m, const uint32_t n){
  for(uint32_t i=0; i<n; i++){
    for(uint32_t k=0; k<m; k++){
      if(k>i){
	if(M[k+i*n] != 0){
	  printf("\nMatrix is not upper triangular; element i:%d j:%d is not zero.\n", i,k);
	  return 0;
	}
      }
    }
  }
  return 0;
}

//Pretty print to user as square array with 5 digits per # with spaces between
//ensuring columns line up...assumes column major storage of elements; prints 
//any sized matrix; must supply m x n dimensions of the matrix
void print_matrix(const double matrix[], const int m, const int n){
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      if((int)(matrix[i+j*m]) >= 100){
	printf("%-9.5f ", matrix[i+j*m]);
      } else if((int)(matrix[i+j*m]) >= 10) {
	printf(" %-9.5f ", matrix[i+j*m]);
      } else {
	printf("  %-9.5f ", matrix[i+j*m]);
      }
    }
    printf("\n");
  }
   printf("\n");
}

