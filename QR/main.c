/* HPC 2012 Project  : Jacqueline Bush, Paul Torres */

/* MAIN FILE: 
Generates n by m matrix and performs Blocked QR factorization/
*/

/* Headers: */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "test.h"
#include "MatrixMatrixMultiply.h"
#include "MatrixTranspose.h"
#include "Utilities.h"


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
  double *D = malloc( n*m*sizeof(double));
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
  //MatrixTranspose(A, m, n, C);
  //MatrixMatrixMultiply( A, n, m, B , m , n , C);  

  for(i=0; i< (m*n) ; i++){
    D[i] = 0;		   
  }
  //MatrixMatrixMultiply( A, n, m, B , m , n , C);  
  
  //printf("=main=A:\n");
  //prettyPrint(A, m, n);

  MatrixTranspose(A, m, n, D);
  
   
  //printf("\nAtranspose:\n");
  //prettyPrint(D, n, m);
   

  //testMatrixTranspose( A, m, n, D);
 
  

}
