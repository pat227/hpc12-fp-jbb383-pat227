/* HPC 2012 Project  : Jacqueline Bush, Paul Torres */

/* BlockedQR2 FILE: 
Performs Blocked QR factorization - version 2
*/

/* Headers: */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "WY.h"
#include "test.h"
#include "Utilities.h"
#include "BlockedQR.h"
#include "BlockedQR2.h"


void BlockedQR2( double *A, int h, int w, double *Q){

/* Block size */
	int b= 16;
	
/* Calculate Number of blocks */
	int wn_bloc = (w+b-1)/b; // Number of Blocks in width (round up)
	int hn_bloc = (h+b -1)/b; // Number of Blocks in height (round up)
	
/* Determine which num. of blocks is smaller */
	int n = wn_bloc; 
	if( wn_bloc > hn_bloc )
		n = hn_bloc;

/* Set Q = m by m identity matrix */
 	for(int i = 0; i<h; i++){
		for(int j=0; j<h; j++){
		Q[i + j*h] = 0;
		if(i == j)
		Q[i +j*h] = 1;
		}
	}	

/* Enter Outer For Loop */
for( int k =0 ; k < n; k++){

  /* Set initial constants for each column */
  int c = 1;
  int c1 = 1; 
  int c2 = 2; 
  int boundary = k;
  
/* Entering While Loop */	
  while(boundary < hn_bloc){
    /* Reset Constants */	
    c = 1;
    /* Compute find the largest c such that for an integer p, c = 2^p < hn_bloc -boundary +1 */ 
    while(c < hn_bloc - boundary + 1){
      c*=2;
    }
    c/=2;
    /* Reset Constants */	
    c1 = 1;
    c2 = 2;
	
  /* Run through Binary Tree - Do in Parallel */
    while( c2 <= c ){
      #pragma omp parallel for shared(A, k, h, w, b, c1, c2, c) 
      for(int i = boundary; i < boundary + c - c1; i += c2){
	//printf("k = %d, BLOCKS called : %d, %d \n", k, i, i+c1);

	/* Call WY Algorithm */
	double  *tempA = malloc(2*b*b*sizeof(double));
	double  *tempR = malloc(2*b*b*sizeof(double));
	double  *tempQ = malloc(2*b*2*b*sizeof(double));
	double *tempQt = malloc(2*b*2*b*sizeof(double));
	Block(tempA, A, h , w, b, i, i + c1, k);	
	WY( tempA, 2*b, b, tempQ, tempQt, tempR);
	UnBlock(A, tempR, h, w, b, i, i + c1, k);
	free(tempA); free(tempR); 

	/* Update Q */
	double *tempQh1 = malloc( b* 2*h* sizeof(double));
	double *tempQh2 = malloc( b* 2*h* sizeof(double));
	BlockQ(Q, tempQh1, h, h, b, i, i + c1 );
	dgemm_simple(tempQh1, h, 2*b, tempQ, 2*b, 2*b, tempQh2);	
	UnBlockQ( tempQh2, Q, h, h, b, i, i+c1);	
	free(tempQ); free(tempQh1); free(tempQh2);
 
	/* Update Blocks along i and i +c1  rows */
	#pragma omp parallel for shared(A, k, h, w, b, c1, i) 	
	for( int j=k+1; j<wn_bloc ; j++){
		double *tempA = malloc(2*b*b*sizeof(double));
		double *tempR = malloc(2*b*b*sizeof(double));
		Block(tempA, A , h, w, b, i, i+c1, j);
		dgemm_simple(tempQt, 2*b, 2*b, tempA, 2*b, b, tempR);
		UnBlock(A, tempR, h, w, b, i, i + c1, j);
		free(tempA); free(tempR);	
	}
	free(tempQt);

      } /* Parallel For Loop Ends */
      
      c1 *= 2; 
      c2 *= 2;
    }
    if((boundary > 0) && (k != boundary)){
      // printf("k = %d, BLOCKS called : %d, %d \n", k, k, boundary);

		double  *tempA = malloc(2*b*b*sizeof(double));
		double  *tempR = malloc(2*b*b*sizeof(double));
		double  *tempQ = malloc(2*b*2*b*sizeof(double));
		double *tempQt = malloc(2*b*2*b*sizeof(double));
		Block(tempA, A, h , w, b, k, boundary, k);	
		WY( tempA, 2*b, b, tempQ, tempQt, tempR);
		UnBlock(A, tempR, h, w, b,k, boundary, k);
		free(tempA); free(tempR); 

		/* Update Q */
		double *tempQh1 = malloc( b* 2*h* sizeof(double));
		double *tempQh2 = malloc( b* 2*h* sizeof(double));
		BlockQ(Q, tempQh1, h, h, b, k, boundary );
		dgemm_simple(tempQh1, h, 2*b, tempQ, 2*b, 2*b, tempQh2);	
		UnBlockQ( tempQh2, Q, h, h, b, k, boundary);	
		free(tempQ); free(tempQh1); free(tempQh2);
 
		/* Update Blocks along boundary and k  rows - Do in Parallel */
		#pragma omp parallel for shared(A, k, h, w, b, boundary) 	
		for( int j=k+1; j<wn_bloc ; j++){
			double *tempA = malloc(2*b*b*sizeof(double));
			double *tempR = malloc(2*b*b*sizeof(double));
			Block(tempA, A , h, w, b, k, boundary, j);
			dgemm_simple(tempQt, 2*b, 2*b, tempA, 2*b, b, tempR);
			UnBlock(A, tempR, h, w, b,k, boundary, j);
			free(tempA); free(tempR);	
		}
		free(tempQt);
    
    } /* End Binary Tree While Loop */
    boundary += c;
    
  } /* Outer While Loop Ends */

/* Update column with only one block in the height */
   if( (k == hn_bloc - 1) && (n == hn_bloc) ){
      //printf("k = %d, BLOCKS called : %d \n", k, k); 

	/* Update Diagonal Block by itself if we are on the last column and there will be no zero rows. */
	double *tempA = malloc( b*b*sizeof(double));
	double *tempR = malloc( b * b* sizeof(double));	
	double *tempQ = malloc( b* b*sizeof(double));
	double *tempQt = malloc( b* b* sizeof(double));
	BlockMatrix(A, tempA, h, w, b, k, k);
	WY( tempA, b, b, tempQ, tempQt, tempR); 
	UnBlockMatrix(A, tempR, h, w, b, k,k);
	free(tempA); free(tempR); 

	/* Update Corresponding Diagonal Block in in Q */
	double *tempQh1 = malloc( b* h* sizeof(double));
	double *tempQh2 = malloc( b* h* sizeof(double));
	BlockQ1(Q, tempQh1, h, h, b, k );
	dgemm_simple(tempQh1, h, b, tempQ, b, b, tempQh2);	
	UnBlockQ1(tempQh2 , Q, h, h, b, k);
	free(tempQ); free(tempQh1); free(tempQh2);
	
	/* Update Blocks along row with Digonal Block - Do in Parallel */
		#pragma omp parallel for shared(A, k, h, w, b, tempQt) 	
		for( int j=k+1; j<wn_bloc ; j++){
		double *tempA = malloc( b*b*sizeof(double));		
		double *tempR = malloc( b*b*sizeof(double));
		BlockMatrix(A, tempA, h, w, b, k, j);
		dgemm_simple(tempQt, b, b, tempA, b, b, tempR);	
		UnBlockMatrix(A, tempR, h,w,b, k, j);
		free(tempA); free(tempR);
		}
   		free(tempQt);
   
    } /* If Loop Ends */

 } /* Outer For Loop Ends */
	

}

