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
/*
  Apply to Q of the A=QR decomposition. Note that Q must be square. I've 
  seen some decompositions called "QR" where Q is not square...these are
  not supported. This function does not store the product matrix, but 
  only computes individual elements at a time and compares against 
  identity matrix.
  m -> # of rows and cols of Q
*/
int IsQbyQtransposeIdentity(const double * const Q, const uint32_t m){
  double temp = 0.0;
  int verbose = 0;
  for(uint32_t i=0; i<m; i++){
    for(uint32_t j=0; j<m; j++){
      for(uint32_t k=0; k<m; k++){
	//QxQt is really computing dot products of all combinations of rows 
	//(since under transpose rows become columns) , and since Q is 
	//orthogonal, they should all = 0 except when dot product is with self
	temp += Q[i+k*m] * Q[j+k*m];
	if(verbose) printf("\ni:%d j:%d k:%d +=  %-9.9f x %-9.9f",
			   i,j,k,Q[i+k*m],Q[j+k*m]);
      }
      if(fabs(temp - identity(i,j)) > EPSILON){
	if(verbose) printf("\nElement i: %d j: %d does not comport with I; computed to be %-9.5f\n", i,j, temp);
	return 0;
      }
      temp = 0.0;
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
int IsQRequalToA(const double * const Q, const double * const R, 
		 const double * const A, const uint32_t Qm, const uint32_t Qn, 
		 const uint32_t Rm, const uint32_t Rn){
  //for each row i and col j in result
  int verbose = 0;
  double temp = 0.0;
  for(uint32_t i=0; i<Qm; i++){
    for(uint32_t j=0; j<Rn; j++){
      for(uint32_t k=0; k<Qn; k++){
	temp += Q[i+k*Qm] * R[k+j*Rm];
	if(verbose) printf("\ni:%d j:%d k:%d +=  %-9.9f x %-9.9f",
			   i,j,k,Q[i+k*Qm],R[k+j*Rm]);
      }
      if(fabs(temp - A[i+j*Qm]) > EPSILON){
	if(verbose) printf("\nA != QR; element i: %d j: %d should be %-9.15f but computed to be %-9.15f\n",
	       i,j,A[i+j*Qm],temp);
	return 0;
      } else {
	temp = 0.0;
      }
    }
  }
  return 1;
}

/*Check constraint that a matrix is upper triangular. Apply to R only in A=QR.*/
int isUpperTriangular(const double * const M, const uint32_t m){
  int verbose = 0;
  for(uint32_t i=0; i<m; i++){
    for(uint32_t j=0; j<i; j++){
      if(fabs(M[i+j*m]) > EPSILON){
	if(verbose) printf("\nMatrix is not upper triangular; element i:%d j:%d is not zero but %f.\n", i,j, M[i+j*m]);
	return 0;
      }
    }
  }
  return 1;
}




//Pretty print to user as square array with 5 digits per # with spaces between
//ensuring columns line up...assumes column major storage of elements; prints 
//any sized matrix; must supply m x n dimensions of the matrix
void print_matrix(const double matrix[], const int m, const int n){
  for(int i = 0; i < m ; i++){
    for(int j = 0; j < n; j++){
      printf("%5.5f ", A[i+j*m]);
    }
    printf("\n");
  }
  printf("\n");
}

NaiveMatrixMultiply(const double * const A, const double * const B, double * C, const int ha, const int wa, const int hb, const int wb){
  for(int i=0; i<mA; i++){
    for(int j=0; j<wB; j++){
      for(int k = 0; k< ; k++){
	C[] = A;
      }
    }
  }  

}
