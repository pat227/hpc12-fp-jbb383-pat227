/* HPC 2012 Project  : Jacqueline Bush, Paul Torres */

/* MAIN FILE: 
Generates n by m matrix and performs Blocked QR factorization/
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include <time.h>
#include "QR.h"
//#include "timing.h"


int main(int argc, char** argv) 
{
/* Check for two arguemnts, m = height of matrix, n = width of matrix.  */ 
  if (argc != 3)
  {
    fprintf(stderr, "Need two arguments, m and n!\n");
    abort(); 
  } 
 
  /* Size of Matrix */  
  int m = atoi(argv[1]);
  int n = atoi(argv[2]);

 /* Initilize matrix A */
  double *A = malloc( m*n*sizeof(double) );
  double *B = malloc( m*n*sizeof(double) );
  double *C = malloc( m*m*sizeof(double) );
  double *Qtranspose = malloc( m*m*sizeof(double) );
 
 /* Fill up matrix A with elements from [0,10)*/
  int i; 
  srand ( time ( NULL ) );	 
  for(i=0; i< (m*n) ; i++){
      A[i] = (double) (rand() %1000)/100;
      B[i] = (double) (rand() %1000)/100;
   }
  for(i=0; i< (m*m) ; i++){
      C[i] = 0;	
      Qtranspose[i] = 0;		   
   }

 // MatrixMatrixMultiply( A, n, m, B , m , n , C);  
  WY( A, n, m, Qtranspose);


}
