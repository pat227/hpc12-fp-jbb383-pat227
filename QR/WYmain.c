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
#include "timing.h"
#include "test.h"
#include "MatrixMatrixMultiply.h"
#include "MatrixTranspose.h"
#include "Utilities.h"
//#include "BlockedQR.h"
#define verbose 1

int main(int argc, char** argv) 
{
/* Check for two arguemnts, m = height of matrix, n = width of matrix k=iterations  */ 
  if (argc != 5)
  {
    fprintf(stderr, "Need three arguments, m, n, iterations, and a (0,1) to indicate if testing is desired.\n");
    abort(); 
  } 
 
  /* Size of Matrix */  
  int m = atoi(argv[1]);
  int n = atoi(argv[2]);
  int iterations = atoi(argv[3]);
  int testing = atoi(argv[4]);
  if(iterations < 1){
    printf("\nIterations must be non-zero positive #.\n");
    abort();
  }

 /* Initilize matrices A, Atest, and Q */
  double * A = (double*) malloc(m * n * sizeof(double));
  double * Atest = (double*) malloc( m * n * sizeof(double));	
  double * Q = (double*) malloc(m * m * sizeof(double));
  double * Qt = (double*) malloc(m * m * sizeof(double));
  double * R = (double *)malloc(m * m * sizeof(double));

 /* Fill up matrix A with elements from [0,10)*/
  int i = 0; 
  srand ( 1 );	 

  timestamp_type time1, time2;
  get_timestamp(&time1); 

  for(int i2 = 0; i2 < iterations; i2++){
    for(i=0; i< (m*n) ; i++){
      A[i] = (double) (rand() %1000)/100;
      Atest[i] = A[i];
    }
    /* BlockedQR replaces A with R, which is why we needed to copy A in order to test code */
    WY(A, m, n, Q, Qt, R);
    if(testing){
      testUpperTriangular(A, m, n);
      printf(" R is Upper Triangular! \n");
      testOrthogonal(Q, Qt, m);
      printf(" Q is Orthogonal! \n");
      IsQRequalToA(Q, A, Atest, m, m, m, n);  	
      printf(" A = QR! QR factorization was sucessful!\n");
    }
  }
  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  double gbs = m * n * 8 * iterations / elapsed / 1e9;

  writetofile("wy_time.txt", m, n, iterations, elapsed);
  writetofile("wy_gbs.txt", m, n, iterations, gbs);
  writetofile2("wy_mbyn.txt", m, n, (elapsed/(double)iterations));
  if(verbose) printf("Time elasped = %f s over %d iterations\n", elapsed, iterations);
  
  free(A);
  free(Atest);
  free(Q);
  free(Qt);
  free(R);

  return 0;
}
