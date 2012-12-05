/* HPC 2012 Project  : Jacqueline Bush, Paul Torres */

/* MAIN FILE: 
Generates n by m matrix and performs Blocked QR factorization/
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include <time.h>
//#include "QR.h"
#include "MatrixMatrixMultiply.h"
//#include "MatrixTranspose.h"

#if 0
//#include "timing.h"
//HAS TO BE STATIC UNTIL .h and .c files re-organized or else get name clash
static void prettyPrint(const double * const A, const int m, const int n){
  printf("m:%d n:%d \n",m,n);
  for(int i = 0; i < m ; i++){
    for(int j = 0; j < n; j++){
      printf("%5.5f ", A[i+j*m]);
    }
    printf("\n");
  }
  printf("\n");
}

#endif

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
  MatrixMatrixMultiply( A, n, m, B , m , n , C);  

  for(i=0; i< (m*n) ; i++){
    D[i] = 0;		   
  }
  //MatrixMatrixMultiply( A, n, m, B , m , n , C);  
  //test2( A, n, m, B, m , n , C);
  printf("=main=A:\n");
  prettyPrint(A, m, n);
  /*

  for(int i = 0; i < m ; i++){
    for(int j = 0; j < n ; j++){
      printf("%5.5f ", A[i+j*m]);
    }
    printf("\n");
  }
  printf("\n");
  */

  //WY( A, n, m, Qtranspose);

  

  printf("\nA again:\n");
  prettyPrint(A, m, n);
  /*
  for(int i = 0; i < m ; i++){
    for(int j = 0; j < n ; j++){
      printf("%5.5f ", A[i+j*m]);
    }
    printf("\n");
  }
  printf("\n");
  */
    
  printf("\nAtranspose:\n");
  prettyPrint(D, n, m);
  /*
  for(int i = 0; i < n ; i++){
    for(int j = 0; j < m ; j++){
      //THIS LINE WAS WRONG AND BAD OUTPUT CAUSED BAD DEBUGGING AND LOTS OF WASTED HOURS; WE NEED FUNCTIONS WITH STANDARDIZED BEHAVIOR with proper .h and .c implementation files!
      printf("%5.5f ", D[i+j*n]);
    }
    printf("\n");
  }
  printf("\n");
  */

  MatrixTranspose(A,n,m,D);
  
  printf("A again:\n");
  for(int i = 0; i < m ; i++){
    for(int j = 0; j < n ; j++){
      printf("%5.5f ", A[i+j*m]);
    }
    printf("\n");
  }
  printf("\n");

  
  printf("Atranspose:\n");
  for(int i = 0; i < m ; i++){
    for(int j = 0; j < n ; j++){
      printf("%5.5f ", D[i+j*m]);
    }
    printf("\n");
  }
  printf("\n");

}
