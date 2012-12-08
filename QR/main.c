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
  double *C = malloc( m*m*sizeof(double) );
 
 
 /* Fill up matrix A with elements from [0,10)*/
  int i; 
  srand ( 1 );	 
  for(i=0; i< (m*n) ; i++){
      A[i] = (double) (rand() %1000)/100;
   }
  for(i=0; i< (m*m) ; i++){
      C[i] = 0;  
   }
  
  WY( A, m, n, C);



   


 

  //for(i=0; i< (m*n) ; i++){
  //  D[i] = 0;		   
 // }
     
  //printf("=main=A:\n");
  //prettyPrint(A, m, n);

  //MatrixTranspose(A, m, n, D);
  
   
  //printf("\nAtranspose:\n");
  //prettyPrint(D, n, m);
   

  //testMatrixTranspose( A, m, n, D);
 
  

}
