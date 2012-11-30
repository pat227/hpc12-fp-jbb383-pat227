#include<stdio.h> //fprintf, printf
#include<stdint.h> //exact types
//#include<stdarg.h> //var args functions
//#include<math.h> //sqrt
#include<stdlib.h> //memory
#include "CUnit/Basic.h" //Basic interface with non-interactive output to stdout.
#include "matrices.h"
extern void init(struct matrix * m, const int w, const int h);
int main(){//int argc, char** argv){
  struct matrix a,b,c,v,h,h2,e1,temp,acopy,q;
  struct matrix * mp;
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
  /*
  a.extractVector(&a,0,0,&v);
  printf("\nV (col vector of A from 0,0):");
  v.prettyPrint(&v);
  norm = norm4(&v,0);
  printf("\nNorm of V: %f", norm);
  e1.zero(&e1);
  if(a.getElement(&a,0,0) < 0){
    e1.setElement(&e1,0,0,1.0);
  } else {
    e1.setElement(&e1,0,0,-1.0);
  }
  e1.scalarMultiply(&e1,norm);
  v.add(&v, &e1, &c);
  v.copyMatrix(&v,&c);
  printf("\nV of Householder:");
  v.prettyPrint(&v);

  v.transpose(&v,&b);
  printf("\nVtranopose:");
  b.prettyPrint(&b);
  printf("\nV x Vtranopose:");
  v.rightMultiply(&v,&b,&c);
  c.prettyPrint(&c);

  printf("\nVtranopose x V:");
  b.rightMultiply(&b, &v, &temp);
  temp.prettyPrint(&temp);

  c.scalarMultiply(&c, 2.0);
  c.scalarMultiply(&c, (1/temp.getElement(&temp, 0, 0)));
  printf("\n 2(v vT / vT v) ");
  c.prettyPrint(&c);
  h.setToIdentity(&h);
  h.subtractFromRightBottomMost(&h, &c);
  printf("\nFirst Householder matrix H(1):");
  h.prettyPrint(&h);

  printf("\nA(1) = H(1) A(0):");
  h.rightMultiply(&h, &a, &temp);
  printf("\nA(1):");
  temp.prettyPrint(&temp);
  a.copyMatrix(&a, &temp);

  //final round
  a.extractVector(&a, 1, 1, &v);
  printf("\nV (col vector of A from 1,1):");
  v.prettyPrint(&v);
  norm = norm4(&v, 0);
  printf("\nNorm of V: %f", norm);
  e1.destroy(&e1);
  e1.init(&e1, 1, 2);
  e1.zero(&e1);
  e1.setElement(&e1, 0, 0, 1.0);
  e1.scalarMultiply(&e1, norm);
  v.add(&v, &e1, &c);
  v.copyMatrix(&v,&c);
  printf("\nV of Householder:");
  v.prettyPrint(&v);

  v.transpose(&v,&b);
  printf("\nVtranopose:");
  b.prettyPrint(&b);
  printf("\nV x Vtranopose:");
  v.rightMultiply(&v,&b,&c);
  c.prettyPrint(&c);

  printf("\nVtranopose x V:");
  b.rightMultiply(&b, &v, &temp);
  temp.prettyPrint(&temp);

  c.scalarMultiply(&c, 2.0);
  c.scalarMultiply(&c, (1/temp.getElement(&temp, 0, 0)));
  printf("\n 2(v vT / vT v) ");
  c.prettyPrint(&c);
  h2.setToIdentity(&h2);
  h2.subtractFromRightBottomMost(&h2, &c);
  printf("\nSecond Householder matrix H(2):");
  h2.prettyPrint(&h2);

  printf("\nA(1) = H(2) A(1):");
  h2.rightMultiply(&h2, &a, &temp);
  printf("\nA(1) (aka R):");
  temp.prettyPrint(&temp);
  a.copyMatrix(&a, &temp);
  
  printf("\nQ (H1 H2) :");
  h.rightMultiply(&h, &h2, &q);
  q.prettyPrint(&q);
  //confirm
  q.rightMultiply(&q,&a,&temp);
  printf("Should equal A:");
  temp.prettyPrint(&temp);
  printf("A again:");
  acopy.prettyPrint(&acopy);
  */
  //looped
  //a.copyMatrix(&acopy);
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
    //mp[i].rightMultiply(&mp[i], &mp[i+1], &h);
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
}
