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
/*
================================================================================
Run tests on known A=QR examples; examples must include where n=m, n<m and n>m
================================================================================
*/
void test_QRisA(){
  //REMEMBER to use COLUMN major storage not Row major
  /*    12 -51   4        6/7  69/175  -58/175       14   21  -14
    A=   6 167 -68   Q =  3/7 -158/175   6/175  R =   0 -175   70
        -4  24 -41       -2/7 -6/35    -33/35         0    0   35
   */   
  //n=m, 3 in this case
  double A[9] = {12,6,-4, -51,167,24,  4,-68,-41};
  double Q[9] = {6.0/7,3.0/7,-2/7.0,   69.0/175,-158/175.0,-6/35.0,
		 -58/175.0,6.0/175,-33/35.0};
  double R[9] = {14,0,0,    21,-175,0,    -14,70,35};
  double wrongR[9] = {14,0,0,    20,-174,0,    -14,70,35};
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
  CU_ASSERT(IsQRequalToA(Q, R, A, 3, 3, 3, 3) == 1);
  CU_ASSERT(IsQRequalToA(Q, wrongR, A, 3, 3, 3, 3) == 0);

  //n<m, 6x3 in this case
  double A2[18] = {4,4,2,1,6,9,   5,8,6,7,1,0,   4,8,6,4,7,2};
  double Q2[18] = {0.322329186, 0.322329186, 0.161164593, 0.080582296, 
		   0.483493778, 0.725240668,
		   0.256776296, 0.513552591, 0.427960493, 0.556348640,
		   -0.171184197,-0.385164443,
		   -0.270918747, 0.148085129, 0.210571917, -0.264284876,
		   0.756475332, -0.467152942};
  double R2[9] = {12.409673646, 0, 0,      6.204836823, 11.683321446, 0,
		  9.992204754, 7.960065161, 4.668319483};
  double wrongQ2[18] = {0.342329186, 0.342329186, 0.161164593, 0.080582296, 
		   0.483493778, 0.725240668,
		   0.256776296, 0.513552591, 0.427960493, 0.556348640,
		   -0.171184197,-0.385164443,
		   -0.270918747, 0.148085129, 0.210571917, -0.264284876,
		   0.756475332, -0.467152942};
  //verbose = 1;
  if(verbose){
    printf("\nOur matrices:");
    printf("\nA:\n");
    print_matrix(A2,6,3);
    printf("\nQ:\n");
    print_matrix(Q2,6,3);
    printf("\nR:\n");
    print_matrix(R2,3,3);
  }
  CU_ASSERT(IsQRequalToA(Q2, R2, A2, 6,3,3,3) == 1);
  CU_ASSERT(IsQRequalToA(wrongQ2, R2, A2, 6,3,3,3) == 0);

  //3x6, n>m case
  double A3[18] = {4,4,2,   5,8,6,   4,8,6,    1,6,9,   7,1,0,    4,7,2};
  double Q3[9] = {0.666666667, 0.666666667, 0.333333333,    
		 -0.630190220, 0.265343251, 0.729693939,      
		 -0.398014876, 0.696526033, -0.597022314};
  double R3[18] = {6.0, 0.0, 0.0,
		   10.666666667, 3.349958540, 0.0,
		   10.0, 3.980148761, 0.398014876,      
		   7.666666667, 7.529114739, -1.592059504,      
		   5.333333333, -4.145988293, -2.089578099,
		   8.0, 0.796029752, 2.089578099};
  //verbose = 1;
  if(verbose){
    printf("\nOur matrices:");
    printf("\nA:\n");
    print_matrix(A3,3,6);
    printf("\nQ:\n");
    print_matrix(Q3,3,3);
    printf("\nR:\n");
    print_matrix(R3,3,6);
  }
  CU_ASSERT(IsQRequalToA(Q3, R3, A3, 3,3,3,6) == 1);
  
}


void test_QbyQtransposeIsIdentity(){
  double Q[9] = {6.0/7,3.0/7,-2/7.0,   69.0/175,-158/175.0,-6/35.0,
		 -58/175.0,6.0/175,-33/35.0};
  /*  THIS CASE FAILS AND BY MY TI-85 as well: QxQt != I, and more to the point, 
      Q is not m x m; however, Qt x Q does = I so perhaps it was bad example
      and not a true QR factorization...but then what is it? */
  //FIXED IT...had been missing 3 more columns; got them from lapack (by way of armadillo++)
  double Q2[36] = {0.322329186, 0.322329186, 0.161164593, 0.080582296, 
		   0.483493778, 0.725240668,
		   0.256776296, 0.513552591, 0.427960493, 0.556348640,
		   -0.171184197,-0.385164443,
		   -0.270918747, 0.148085129, 0.210571917, -0.264284876,
		   0.756475332, -0.467152942,
		   -0.17063541, -0.632751924, 0.145160343, 0.678761105, 
		   0.292902314, 0.054117043,
		   -0.182799363, 0.358320825, -0.806791034, 0.386289771, 
		   0.184148016, -0.0644090721,
		   -0.833211159, 0.285778962, 0.273042051, 0.0648055884,
		   -0.212067327, 0.316804785  };
    double Q3[9] = {0.666666667, 0.666666667, 0.333333333,    
		 -0.630190220, 0.265343251, 0.729693939,      
		 -0.398014876, 0.696526033, -0.597022314};
    double wrongQ3[9] = {0.666666667, 0.333333333, 0.666666667,
		 -0.630190220, 0.265343251, 0.729693939,      
		 -0.398014876, 0.696526033, -0.597022314};
    int verbose = 0;
  if(verbose){
    printf("\nOur matrices:");
    printf("\nQ:\n");
    print_matrix(Q,3,3);
    printf("\nQ2:\n");
    print_matrix(Q2,6,6);
    printf("\nQ3:\n");
    print_matrix(Q3,3,3);
  }  
  CU_ASSERT(IsQbyQtransposeIdentity(Q,3) == 1);
  CU_ASSERT(IsQbyQtransposeIdentity(Q2,6) == 1);
  CU_ASSERT(IsQbyQtransposeIdentity(Q3,3) == 1);
  CU_ASSERT(IsQbyQtransposeIdentity(wrongQ3,3) == 0);
}

void test_RisUpperTriangular(){
  double R[9] = {14,0,0,    21,-175,0,    -14,70,35};
  double R2[9] = {12.409673646, 0, 0,      6.204836823, 11.683321446, 0,
		  9.992204754, 7.960065161, 4.668319483};
  double R3[18] = {6.0, 0.0, 0.0,
		   10.666666667, 3.349958540, 0.0,
		   10.0, 3.980148761, 0.398014876,      
		   7.666666667, 7.529114739, -1.592059504,      
		   5.333333333, -4.145988293, -2.089578099,
		   8.0, 0.796029752, 2.089578099};
  double wrongR3[18] = {6.0, 0.0, 0.0,
		   10.666666667, 3.349958540, 0.0001,
		   10.0, 3.980148761, 0.398014876,      
		   7.666666667, 7.529114739, -1.592059504,      
		   5.333333333, -4.145988293, -2.089578099,
		   8.0, 0.796029752, 2.089578099};
  CU_ASSERT(isUpperTriangular(R,3) == 1);
  CU_ASSERT(isUpperTriangular(R2,3) == 1);
  CU_ASSERT(isUpperTriangular(R3,3) == 1);
  CU_ASSERT(isUpperTriangular(wrongR3,3) == 0);
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
   if (NULL == CU_add_test(pSuite, "id()", test_id) || 
       NULL == CU_add_test(pSuite,"A=QR",test_QRisA) || 
       NULL == CU_add_test(pSuite, "QxQtranspose=I", test_QbyQtransposeIsIdentity) ||
       NULL == CU_add_test(pSuite, "R is upper Triangular", test_RisUpperTriangular)){
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
