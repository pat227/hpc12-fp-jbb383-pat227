/* HPC 2012 Project  : Jacqueline Bush, Paul Torres */

/* MAIN FILE: 
Generates n by m matrix and performs Blocked QR factorization/
*/

/* Headers: */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "WY.h"
#include "test.h"
#include "MatrixMatrixMultiply.h"
#include "MatrixTranspose.h"
#include "Utilities.h"
#include "BlockedQR.h"


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
  double *B = malloc( m*n*sizeof(double));	
  double *C = malloc( m*m*sizeof(double));	

  double *Q = malloc( m*m*sizeof(double) );
  double *R = malloc(m *n *sizeof(double));
  
 
 
 /* Fill up matrix A with elements from [0,10)*/
  int i; 
  srand ( 1 );	 
  for(i=0; i< (m*n) ; i++){
      A[i] = (double) (rand() %1000)/100;
      B[i] = (double) (rand() %1000)/100;
   }
  
 MatrixMatrixMultiply(A, m, n, B, n, m, C);
 
  BlockedQR( A, m, n, Q, R);

  

}
