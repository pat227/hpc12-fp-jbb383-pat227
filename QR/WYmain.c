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
#include "BlockedQR.h"
#define verbose 1

int main(int argc, char** argv) 
{
/* Check for two arguemnts, m = height of matrix, n = width of matrix k=iterations  */ 
  if (argc != 4)
  {
    fprintf(stderr, "Need three arguments, m, n, iterations.\n");
    abort(); 
  } 
 
  /* Size of Matrix */  
  int m = atoi(argv[1]);
  int n = atoi(argv[2]);
  int iterations = atoi(argv[3]);

 /* Initilize matrices A, Atest, and Q */
  double * A = malloc(m * n * sizeof(double));
  double * Atest = malloc( m * n * sizeof(double));	
  double * Q = malloc(m * m * sizeof(double));
  double * Qt = malloc(m * m * sizeof(double));
  double * R = malloc(m * m * sizeof(double));

 /* Fill up matrix A with elements from [0,10)*/
  int i = 0; 
  srand ( 1 );	 

  timestamp_type time1, time2;
  get_timestamp(&time1); 

  for(int i2 = 0; i2 < iterations; i2++){
    for(i=0; i< (m*n) ; i++){
      A[i] = (double) (rand() %1000)/100;
      //Atest[i] = A[i];
    }
    WY(A, m, n, Q, Qt, R);
  }
  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  double gbs = m * n * 8 * iterations / elapsed / 1e9;

  writetofile("wy_time.txt", m, n, iterations, elapsed);
  writetofile("wy_gbs.txt", m, n, iterations, gbs);
  
  if(verbose) printf("Time elasped = %f s over %d iterations\n", elapsed, iterations);

  free(A);
  free(Atest);
  free(Q);
  free(Qt);
  free(R);

  return 0;
}
