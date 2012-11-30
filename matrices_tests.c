#include<stdio.h> //fprintf, printf
//#include<stdint.h> //exact types
//#include<stdarg.h> //var args functions
//#include<math.h> //sqrt
//#include<stdlib.h> //memory
#include "CUnit/Basic.h" //Basic interface with non-interactive output to stdout.
#include "matrices.h"
extern void init(struct matrix * m, const int w, const int h);
void test_BasicFunctions(){
  struct matrix a;
  struct matrix b;
  struct matrix c;
  struct matrix * mp;
  printf("\nInitializing...\n");
  init(&a, 3,3);
  init(&b, 3,3);
  //alternative style...use the pointer
  mp = &c;
  init(mp, 3,3);

  zero(&a);
  for(int i = 0; i < a.height; i++){
    for(int j = 0; j < a.width; j++){
      CU_ASSERT(getElement(&a,i,j) == 0.0);
    }
  }

  setToIdentity(&b);
  for(int i = 0; i < a.height; i++){
    for(int j = 0; j < a.width; j++){
      if(i==j){
	CU_ASSERT(getElement(&b,i,j) == 1.0);
      } else {
	CU_ASSERT(getElement(&b,i,j) == 0.0);
      }
    }
  }
  zero(&c);
  
  setElement(&a, 0,0,12.0); 
  CU_ASSERT(getElement(&a,0,0) == 12.0);
  setElement(&a, 1,0,6.0);
  setElement(&a, 2,0,-4.0);
  setElement(&a, 0,1,-51.0);
  setElement(&a, 1,1,167.0);
  CU_ASSERT(getElement(&a,1,1) == 167.0);
  setElement(&a, 2,1,24.0);
  setElement(&a, 0,2,4.0);
  setElement(&a, 1,2,-68.0);
  CU_ASSERT(getElement(&a,1,2) == -68.0);
  setElement(&a, 2,2,-41.0);
  printf("A:");
  prettyPrint(&a);
  printf("B:");
  prettyPrint(&b);

  copyMatrix(&c,&a);  
  CU_ASSERT(getElement(&c,1,1) == getElement(&a,1,1));
  printf("C (copy of A):");
  prettyPrint(&c);

}
void test_moreFunctions(){
  struct matrix a,b,c,v,h,h2,e1,temp,acopy,q;
  //struct matrix * mp;
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
  //  b.setToIdentity(&b);
  //a.zero(&a);
  setElement(&a, 0,0,12.0); 
  setElement(&a, 1,0,6.0);
  setElement(&a, 2,0,-4.0);
  setElement(&a, 0,1,-51.0);
  setElement(&a, 1,1,167.0);
  setElement(&a, 2,1,24.0);
  setElement(&a, 0,2,4.0);
  setElement(&a, 1,2,-68.0);
  setElement(&a, 2,2,-41.0);
  copyMatrix(&acopy, &a);
  printf("A:");
  prettyPrint(&a);
  printf("Should match:\n");
  for(int i = 0; i < a.height; i++){
    for(int j = 0; j < a.width; j++){
      printf(" %f ", getElement(&a, i,j));
    }
    printf("\n");
  }
  
  printf("B:");
  setToIdentity(&b);
  prettyPrint(&b);

  matrixMultiply(&a,&b,&c);
  printf("C (=AB):");
  prettyPrint(&c);
  printf("C (Asquared):");
  matrixMultiply(&a,&b,&c);
  prettyPrint(&c);
  zero(&c);

  scalarMultiply(&a,2.0);
  printf("A (x2):");
  prettyPrint(&a);
  scalarMultiply(&a,0.5);
  printf("A (/2):");
  prettyPrint(&a);
  add(&a, &b, &c);

  printf("B again:");
  prettyPrint(&b);
  printf("C (A+B):");
  prettyPrint(&c);
  scalarMultiply(&b,-1.0);
  add(&a, &b, &c);
  printf("C (A-B):");
  prettyPrint(&c);

  printf("A again:");
  prettyPrint(&a);
  transpose(&a,&c);
  printf("C (tranpose of A):");
  prettyPrint(&c);

  extractVector(&a,0,0,&v);
  printf("V (first col vector of A):");
  prettyPrint(&v);
  printf("\nNorm of V: %f", normOfVector(&v,0));
  extractVector(&a,0,1,&v);
  printf("\nV (2nd col vector of A):");
  prettyPrint(&v);
  printf("\nNorm of V: %f", normOfVector(&v,0));
  extractVector(&a,0,2,&v);
  printf("\nV (third col vector of A):");
  prettyPrint(&v);
  printf("\nNorm of V: %f", normOfVector(&v,0));   
  printf("\nNorm of V(starting at 2nd row & col): %f", normOfColVector(&a,1,1));

  printf("V (vector of A from 1,0):");
  extractVector(&a,1,0,&v);
  prettyPrint(&v);
  printf("\nNorm of V: %f", normOfVector(&v,0));
  printf("\nV (col vec of A from 1,1):");
  extractVector(&a,1,1,&v);
  prettyPrint(&v);
  printf("\nNorm of V: %f", normOfVector(&v,0));
  extractVector(&a,1,2,&v);
  printf("\nV (third col vector of A from 1,3):");
  prettyPrint(&v);
  printf("\nNorm of V: %f", normOfVector(&v,0));   
  printf("\nNorm of V(starting at 2nd row & col): %f", normOfColVector(&a,1,1));

}

int main(){//int argc, char** argv){
  CU_pSuite pSuite = NULL;
  
  /* initialize the CUnit test registry */
  if (CUE_SUCCESS != CU_initialize_registry())
    return CU_get_error();
  
  /* add a suite to the registry */
  pSuite = CU_add_suite("Suite_1", NULL,NULL);//init_suite1, clean_suite1);
  if (NULL == pSuite) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  
  /* add the tests to the suite */
  if (NULL == CU_add_test(pSuite, "BasicFunctions", test_BasicFunctions) || 
      NULL == CU_add_test(pSuite, "MoreFunctions", test_moreFunctions)){
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  //CU_console_run_tests();
  CU_cleanup_registry();
  return CU_get_error();
}
