//#include <CUnit/CUnit.h>	//ASSERT macros for use in test cases, and includes other framework headers.
//#include <CUnit/CUError.h>	//Error handing functions and data types. Included automatically by CUnit.h.
//#include <CUnit/TestDB.h>	//Data type definitions and manipulation functions for the test registry, suites, and tests. Included automatically by CUnit.h.
//#include <CUnit/TestRun.h>	//Data type definitions and functions for running tests and retrieving results. Included automatically by CUnit.h.
//#include <CUnit/Automated.h>	//Automated interface with xml output.
#include "CUnit/Basic.h"	//Basic interface with non-interactive output to stdout.
//#include <CUnit/Console.h>	//Interactive console interface.
//#include <CUnit/CUCurses.h>	//Interactive console interface (*nix).

#include<stdio.h> //fprintf, printf
#include<stdint.h> //exact types
#include"proj-LibCorrectness.h" //the functions we wish to test have headers here
/*We use column major access and assume args are stored in column major order in accordance with all the other work we did this semester*/


void test_id(void){
  CU_ASSERT(identity(0,2) == 0);
  CU_ASSERT(identity(0,-2) == 0);
  CU_ASSERT(identity(2,2) == 1);
  CU_ASSERT(identity(1,1) == 1);
  CU_ASSERT(identity(0,0) == 1);
}

void test_QRisA(){
  //REMEMBER to use COLUMN major storage not Row major
  /*    12 -51   4        6/7  69/175  -58/175       14   21  -14
    A=   6 167 -68   Q =  3/7 -158/175   6/175  R =   0 -175   70
        -4  24 -41       -2/7 -6/35    -33/35         0    0   35
   */   
  double A[9] = {12,6,-4, -51,167,24,  4,-68,-41};
  double Q[9] = {6.0/7,3.0/7,-2/7.0,   69.0/175,-158/175.0,-6/35.0,  -58/175.0,6.0/175,-33/35.0};
  double R[9] = {14,0,0,    21,-175,0,    -14,70,35};
  int verbose = 0;
  if(verbose){
    printf("\nOur matrices:");
    printf("\nA:\n");
    print_matrix(A,3,3);
    printf("\nQ:\n");
    print_matrix(Q,3,3);
    printf("\nR:\n");
    print_matrix(R,3,3);
  }
  CU_ASSERT(IsQRequalToA(Q, R, A, 3) == 1);
  //do some more here...
  /*
  double A2[9] = {12,6,-4, -51,167,24,  4,-68,-41};
  double Q2[9] = {6.0/7,3.0/7,-2/7.0,   69.0/175,-158/175.0,-6/35.0,  -58/175.0,6.0/175,-33/35.0};
  double R2[9] = {14,0,0,    21,-175,0,    -14,70,35};
  int verbose = 0;
  if(verbose){
    printf("\nOur matrices:");
    printf("\nA:\n");
    print_matrix(A2,3,3);
    printf("\nQ:\n");
    print_matrix(Q2,3,3);
    printf("\nR:\n");
    print_matrix(R2,3,3);
  }
  CU_ASSERT(IsQRequalToA(Q2, R2, A2, 3) == 1);
  */
}

int main()
{
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
   if (NULL == CU_add_test(pSuite, "id()", test_id) || NULL == CU_add_test(pSuite,"A=QR",test_QRisA)){
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
