/* HPC 2012 Project  : Jacqueline Bush, Paul Torres */

/* MAIN FILE: 
Generates n by m matrix and performs Blocked QR factorization/
*/

/* Headers: */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "timing.h"
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

 /* Initilize matrices A, Atest, and Q */
  double *A = malloc( m*n*sizeof(double) );
  double *Atest = malloc(m*n*sizeof(double));	
  double *Q = malloc( m*m*sizeof(double) );

 /* Fill up matrix A with elements from [0,10)*/
  int i; 
  srand ( 1 );	 
  for(i=0; i< (m*n) ; i++){
      A[i] = (double) (rand() %1000)/100;
      Atest[i] = A[i];	
   }
 
  /* Start Timing */
  timestamp_type time1, time2;
  get_timestamp(&time1); 

 /* BlockedQR replaces A with R, which is why we needed to copy A in order to test code */
  BlockedQR( A, m, n, Q);

 /* End Timing */
  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  printf("Blocked QR takes %f s\n", elapsed);
 
 /* Test Code */
 printf(" Testing Code ..... \n");
 
 double *Qt = malloc( m* m* sizeof(double)); 
 MatrixTranspose(Q, m, m, Qt);	

 testUpperTriangular(A, m, n);
 printf(" R is Upper Triangular! \n");
 testOrthogonal(Q, Qt, m);
 printf(" Q is Orthogonal! \n");
 IsQRequalToA(Q, A, Atest, m, m, m, n);  	
 printf(" A = QR! QR factorization was sucessful!\n");

 /* Clean Up */
  free(A); free(Atest); free(Q); free(Qt);

}
