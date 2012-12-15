#include<stdio.h> //fprintf, printf
#include<stdint.h> //exact types
#include<stdlib.h> //memory
#include "matrices.h"
#include "proj-LibCorrectness.h"
#include "timing.h"
#define verbose 0

extern void init(struct matrix * m, const int w, const int h);
int main(int argc, char** argv){
  if(argc != 5){
    printf("\n================================================================================");
    printf("\nGeneralized householder & example arguments:                                    ");
    printf("\nm -> desired height of matrix to decompose into QR                              ");
    printf("\nn -> desired width of matrix to decompose into QR                               ");
    printf("\niterations -> desied # of iterations (for performance timing purposes)          ");
    printf("\ntesting -> (0,1) to indicate if testing is desired.");
    printf("\nThis program will compute an example 3x3 QR decomposition by Householder        ");
    printf("\nreflectors and then compute a QR decomposition for a matrix of random elements  ");
    printf("\nof the size specified at the command line by the 2 preceding arguments.         ");
    printf("\nHPC Fall'12 GSAS NYU                                                            ");
    printf("\n================================================================================\n");
    return 1;
  }

  struct matrix a,b,c,v,h,h2,e1,temp,acopy,q;
  struct matrix * mp = NULL;
  int number = 0;
  int j = 0;
  double norm = 0;
  int m = atoi(argv[1]);
  int n = atoi(argv[2]);
  if(m == 0 || n == 0){
    printf("\nCould not parse the desired matrix dimensions m, n");
    abort();
  } else if(m < 0 || n < 0){
    printf("\nMatrix dimensions m, n must be greater than zero.");
    abort();
  }
  int iterations = atoi(argv[3]);
  int testing = atoi(argv[4]);
  if(iterations < 1){
    printf("\nIterations must be a non-zero positive number.\n");
  }

  timestamp_type time1, time2;
  get_timestamp(&time1);

  for(int iter = 0; iter < iterations; iter++){
    init(&b, m,n);
    init(&c, m,n);
    init(&v, m,1);
    init(&h, m,m);
    init(&h2, m,m);
    init(&e1, m,1);
    init(&temp, m,n);
    init(&q, m,m);
    //============================================================================
    //generalized version using supplied args
    //============================================================================
    mp = NULL;
    number = m;
    if(m > n){
      number = n;//(n-1);
    }
    mp = malloc(sizeof(struct matrix)*(m+1));
    if(mp == NULL){
      printf("Couldn't allocate matrices array.");
      abort();
    }
    if(verbose){
      printf("\n===============================================================================");
      printf("\nPerforming QR decomposition upon random matrix of size m x n...");
      printf("\n===============================================================================");
    }
    for(int i = 0; i < (m+1); i++){
      init(&mp[i], m, m);
      setToIdentity(&mp[i]);
    }
    //reset j and a
    j = 0;
    init(&a, n, m);
    init(&acopy, n, m);
    fillWithRandomElements(&a, 10, 0);
    if(verbose){
      printf("A:");
      prettyPrint(&a);
    }
    copyMatrix(&acopy, &a);
    
    //the loop --AGAIN -- that computes all the reflectors and QR
    while(j+1 <= a.width && j+1 <= a.height){
      extractVector(&a,j,j,&v);
      
      if(verbose){
	printf("\nV (col vector of A from j,j):");
	prettyPrint(&v);
      }
      
      norm = normOfVector(&v,0);
      
      if(verbose) printf("\nNorm of V: %f", norm);
      
      init(&e1, 1, a.height-j);
      zero(&e1);
      if(getElement(&a,j,j) < 0){
	setElement(&e1,0,0,-1.0);
      } else {
	setElement(&e1,0,0,1.0);
      }
      scalarMultiply(&e1, norm);
      add(&v, &e1, &c);
      //use swap instead
      //copyMatrix(&v, &c);
      swapMatrix(&v, &c);
      
      if(verbose) {
	printf("\nV of Householder:");
	prettyPrint(&v);
      }

      transpose(&v,&b);
      if(verbose){
	printf("\nVtranopose:");
	prettyPrint(&b);
	printf("\nV x Vtranopose:");
      }
      matrixMultiply(&v,&b,&c);
      
      if(verbose){
	prettyPrint(&c);
	printf("\nVtranopose x V:");
      }
      matrixMultiply(&b, &v, &temp);
      if(verbose) prettyPrint(&temp);
      
      scalarMultiply(&c, 2.0);
      scalarMultiply(&c, (1 / getElement(&temp, 0, 0)));
      if(verbose){
	printf("\n 2(v vT / vT v) ");
	prettyPrint(&c);
      }
      //setToIdentity(&(mp[j]));
      subtractFromRightBottomMost(&mp[j], &c);
      if(verbose){
	printf("\n%d Householder matrix H(%d):",j+1, j+1);
	prettyPrint(&mp[j]);
      }
      
      if(verbose) printf("\nA(%d) = H(%d) A(%d):",j+1,j+1,j);
      matrixMultiply(&mp[j], &a, &temp);
      if(verbose){
	printf("\nA(%d):",j+1);
	prettyPrint(&temp);
      }
      //use swap instead
      //copyMatrix(&a, &temp);    
      swapMatrix(&a, &temp);
      j++;
    }  
    init(&h, m, m);
    setToIdentity(&h);
    
    //printf("j bound: %d ", j);
    for(int i = 0; i < j; i+=2){
      matrixMultiply(&h, &mp[i], &h2);
      //use swap instead
      //copyMatrix(&h, &h2);
      swapMatrix(&h,&h2);
      matrixMultiply(&h, &mp[i+1], &h2);
      //use swap instead
      //copyMatrix(&h, &h2);
      swapMatrix(&h,&h2);
    }
    copyMatrix(&q, &h);
    //confirm -- use functions not on screen output
    if(verbose){
      printf("\nQ:");
      prettyPrint(&q);
      matrixMultiply(&q, &a, &temp);
      printf("Should equal A:");
      prettyPrint(&temp);
      printf("A again:");
      prettyPrint(&acopy);
    }
    if(verbose) printf("\nChecking that QR=A, QQtranspose = I, and that R is upper triangular...");
    //copyMatrix(&q, &h);
    if(testing){
      if(IsQRequalToA(q.elements, a.elements, acopy.elements, m, m, m, n)){
	printf("\nQR = A checks...");
      } else {
	printf("\nQR = A DOES NOT check...");
      }
      if(IsQbyQtransposeIdentity(h.elements, m)){
	printf("\nQQtranspose = I checks...");
      } else {
	printf("\nQQtranspose = I DOES NOT check...");
      }
      if(isUpperTriangular(a.elements, m)){
	printf("\nR is upper triangular checks...");
      } else {
	printf("\nR is upper triangular DOES NOT  check...");
      }
    }
    
  }

  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2);
  printf("Total Elapsed Time: %f", elapsed);
  
  printf("\n");

  //============cleanup=============
  free(mp);
  return 0;
}
