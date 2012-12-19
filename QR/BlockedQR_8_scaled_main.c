/* HPC 2012 Project  : Jacqueline Bush, Paul Torres */

/* MAIN FILE: 
Generates n by m matrix and performs Blocked QR factorization/
*/

/* Headers: */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
//#include "WY.h"
#include "timing.h"
#include "test.h"
#include "MatrixMatrixMultiply.h"
#include "MatrixTranspose.h"
#include "Utilities.h"
#include "BlockedQR.h"
#define verbose 1

int main(int argc, char** argv) 
{
/* Check for arguemnts, m = height of matrix, n = width of matrix k=iterations  */ 
  if (argc != 5)
  {
    fprintf(stderr, "Need four arguments, m, n, iterations, and a (0,1) to indicate if testing is desired.\n");
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
  double * A = malloc(m * n * sizeof(double));
  double * Atest = malloc( m * n * sizeof(double));	
  double * Q = malloc(m * m * sizeof(double));

 /* Fill up matrix A with elements from [0,10)*/
  int i = 0; 
  srand ( time(NULL) );	 

  timestamp_type time1, time2;
  get_timestamp(&time1); 

  for(int i2 = 0; i2 < iterations; i2++){
    for(i=0; i< (m*n) ; i++){
      A[i] = (double) (rand() %1000)/100;
      Atest[i] = A[i];
    }
    /* BlockedQR replaces A with R, which is why we needed to copy A in order to test code */
    BlockedQR(A, m, n, Q);

    /* Test Code */
    if(testing){
      printf(" Testing Code ..... \n");
      
      double *Qt = malloc( m* m* sizeof(double)); 
      MatrixTranspose(Q, m, m, Qt);	
      
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
  //double gbs = m * n * 8 * iterations / elapsed / 1e9;
  double kflops = m * n * n / elapsed / 1000; //fabs(m * n * n / 2 - n * n * n / 3);
  char * ompthreads = getenv("NUM_OMP_THREADS");
  printf("\n ran with OMP %s ", ompthreads);
  //just use output...bug someplace...and not many points to copy paste
  //writetofile3("blockedQR_8_scaled_time.txt", elapsed, ompthreads);
  //writetofile3("blockedQR_8_scaled_gfs.txt", kflops, ompthreads);
  
  if(verbose) printf("QR_8_scaled: Time = %f numthreads=%s, m:%d, n:%d kflops:%f\n", elapsed, ompthreads, m,n, kflops);

  free(A);
  free(Atest);
  free(Q);

  return 0;
}
