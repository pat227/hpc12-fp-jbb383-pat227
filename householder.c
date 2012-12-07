#include<stdio.h> //fprintf, printf
#include<stdint.h> //exact types
#include<stdlib.h> //memory
#include "matrices.h"
#include "proj-LibCorrectness.h"
#define verbose 1
extern void init(struct matrix * m, const int w, const int h);
int main(int argc, char** argv){
  if(argc != 3){
    printf("\n================================================================================");
    printf("\nGeneralized householder & example arguments:                                    ");
    printf("\nm -> desired height of matrix to decompose into QR                              ");
    printf("\nn -> desired width of matrix to decompose into QR                               ");
    printf("\nThis program will compute an example 3x3 QR decomposition by Householder        ");
    printf("\nreflectors and then compute a QR decomposition for a matrix of random elements  ");
    printf("\nof the size specified at the command line by the 2 preceding arguments.         ");
    printf("\nHPC Fall'12 GSAS NYU                                                            ");
    printf("\n================================================================================\n");
    return 1;
  }

  struct matrix a,b,c,v,h,h2,e1,temp,acopy,q;
  struct matrix * mp;
  int m = atoi(argv[1]);
  int n = atoi(argv[2]);
  if(m == 0 || n == 0){
    printf("\nCould not parse the desired matrix dimensions m, n");
    abort();
  } else if(m < 0 || n < 0){
    printf("\nMatrix dimensions m, n must be greater than zero.");
    abort();
  }
  init(&a, 3,3);
  init(&b, 3,3);
  init(&c, 3,3);
  //can make it anything actually, since it will reshape itself when I extract a column vector later
  init(&v, 3,3);
  init(&h, 3,3);
  init(&h2, 3,3);
  init(&e1, 1,3);
  init(&temp,3,3);
  init(&acopy, 3,3);
  init(&q ,3,3);

  setToIdentity(&b);
  zero(&a);
  setElement(&a, 0,0,12.0); 
  setElement(&a, 1,0,6.0);
  setElement(&a, 2,0,-4.0);
  setElement(&a, 0,1,-51.0);
  setElement(&a, 1,1,167.0);
  setElement(&a, 2,1,24.0);
  setElement(&a, 0,2,4.0);
  setElement(&a, 1,2,-68.0);
  setElement(&a, 2,2,-41.0);
  printf("A:");
  prettyPrint(&a);
  copyMatrix(&acopy, &a);
  int number = 0;
  if(a.height > a.width){
    number = a.height-1;
  } else {
    number = a.width-1;
  }
  mp = malloc(sizeof(struct matrix)*number);
  if(mp == NULL){
    printf("Couldn't allocate matrices array.");
    abort();
  }
  for(int i = 0; i < number; i++){
    init(&mp[i], 3, 3);
  }
  //============================================================================
  //looped
  //============================================================================
  printf("=Looped=");
  double norm = 0.0;
  int j = 0;
  while(j+1 < a.width || j+2 < a.height){
    extractVector(&a,j,j,&v);
    printf("\nV (col vector of A from j,j):");
    prettyPrint(&v);
    norm = normOfVector(&v,0);
    printf("\nNorm of V: %f", norm);
    init(&e1, 1, a.height-j);
    zero(&e1);
    if(getElement(&a,j,j) < 0){
      setElement(&e1,0,0,1.0);
    } else {
      setElement(&e1,0,0,-1.0);
    }  
    scalarMultiply(&e1, norm);
    add(&v, &e1, &c);
    copyMatrix(&v, &c);
    printf("\nV of Householder:");
    prettyPrint(&v);

    transpose(&v,&b);
    printf("\nVtranopose:");
    prettyPrint(&b);
    printf("\nV x Vtranopose:");
    matrixMultiply(&v,&b,&c);
    prettyPrint(&c);
    
    printf("\nVtranopose x V:");
    matrixMultiply(&b, &v, &temp);
    prettyPrint(&temp);
    
    scalarMultiply(&c, 2.0);
    scalarMultiply(&c, (1 / getElement(&temp, 0, 0)));
    printf("\n 2(v vT / vT v) ");
    prettyPrint(&c);
    setToIdentity(&(mp[j]));
    subtractFromRightBottomMost(&mp[j], &c);
    printf("\n%d Householder matrix H(%d):",j+1, j+1);
    prettyPrint(&mp[j]);

    printf("\nA(%d) = H(%d) A(%d):",j+1,j+1,j);
    matrixMultiply(&mp[j], &a, &temp);
    printf("\nA(%d):",j+1);
    prettyPrint(&temp);
    copyMatrix(&a, &temp);    
    j++;
  }  
  setToIdentity(&h);
  for(int i = 0; i < j; i+=2){
    matrixMultiply(&h, &mp[i], &h2);
    copyMatrix(&h, &h2);
    matrixMultiply(&h, &mp[i+1], &h2);
    copyMatrix(&h, &h2);
  }
  printf("\nQ:");
  copyMatrix(&q, &h);
  prettyPrint(&q);
  //confirm
  matrixMultiply(&q, &a, &temp);
  printf("Should equal A:");
  prettyPrint(&temp);
  printf("A again:");
  prettyPrint(&acopy);  

  //============================================================================
  //generalized version using supplied args
  //============================================================================
  free(mp);
  mp = NULL;
  number = 0;
  if(m > n){
    number = (n-1);
  } else {
    number = (m-1);
  }
  mp = malloc(sizeof(struct matrix)*number);
  if(mp == NULL){
    printf("Couldn't allocate matrices array.");
    abort();
  }
  printf("\n===============================================================================");
  printf("\nPerforming QR decomposition upon random matrix of size m x n...");
  printf("\n===============================================================================");
  for(int i = 0; i < number; i++){
    init(&mp[i], m, m);
  }
  //reset j and a
  j = 0;
  init(&a, n, m);
  fillWithRandomElements(&a, 10, 0);
  if(verbose){
    printf("A:");
    prettyPrint(&a);
  }
  copyMatrix(&acopy, &a);
  //the loop --AGAIN -- that computes all the reflectors and QR
  while(j+1 < a.width || j+2 < a.height){
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
      setElement(&e1,0,0,1.0);
    } else {
      setElement(&e1,0,0,-1.0);
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
    setToIdentity(&(mp[j]));
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
  //confirm -- use functions not on screen output
  if(verbose){
    printf("\nQ:");
    copyMatrix(&q, &h);
    prettyPrint(&q);
    matrixMultiply(&q, &a, &temp);
    printf("Should equal A:");
    prettyPrint(&temp);
    printf("A again:");
    prettyPrint(&acopy); 
  }
  printf("\nChecking that QR=A, QQtranspose = I, and that R is upper triangular...");
  //copyMatrix(&q, &h);
  if(IsQRequalToA(h.elements, a.elements, acopy.elements, m, m, n, n)){
    printf("\nQR = A checks...");
  }
  if(IsQbyQtransposeIdentity(h.elements, m)){
    printf("\nQQtranspose = I checks...");
  }
  if(isUpperTriangular(a.elements, m)){
    printf("\nR is upper triangular checks...");
    }
  printf("\n");
  return 0;
}
