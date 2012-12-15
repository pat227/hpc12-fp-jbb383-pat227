
/* HPC 2012 Project  : Jacqueline Bush, Paul Torres */

/* BlockedQR FILE: 
Performs Blocked QR factorization/
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
	int b= 16 ;

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


/* Enter Loop */
for(int k =0; k< n; k++){	
	/* Update Diagonal Block */
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

	/* Update Blocks below the Diagonal (will also update diagonal) */
	for( int i=k+1; i<hn_bloc ; i++){
		double  *tempA = malloc(2*b*b*sizeof(double));
		double  *tempR = malloc(2*b*b*sizeof(double));
		double  *tempQ = malloc(2*b*2*b*sizeof(double));
		double *tempQt = malloc(2*b*2*b*sizeof(double));
		Block(tempA, A, h , w, b, k, i, k);	
		WY( tempA, 2*b, b, tempQ, tempQt, tempR);
		UnBlock(A, tempR, h, w, b,k, i, k);
		free(tempA); free(tempR); 

		/* Update Q */
		double *tempQh1 = malloc( b* 2*h* sizeof(double));
		double *tempQh2 = malloc( b* 2*h* sizeof(double));
		BlockQ(Q, tempQh1, h, h, b, k, i );
		dgemm_simple(tempQh1, h, 2*b, tempQ, 2*b, 2*b, tempQh2);	
		UnBlockQ( tempQh2, Q, h, h, b, k, i);	
		free(tempQ); free(tempQh1); free(tempQh2);
 
		/* Update Blocks along i and k  rows */
		#pragma omp parallel for shared(A, k, h, w, b, i) 	
		for( int j=k+1; j<wn_bloc ; j++){
			double *tempA = malloc(2*b*b*sizeof(double));
			double *tempR = malloc(2*b*b*sizeof(double));
			Block(tempA, A , h, w, b, k, i, j);
			dgemm_simple(tempQt, 2*b, 2*b, tempA, 2*b, b, tempR);
			UnBlock(A, tempR, h, w, b,k, i, j);
			free(tempA); free(tempR);	
		}
		free(tempQt);

	}

	
}

}






 









