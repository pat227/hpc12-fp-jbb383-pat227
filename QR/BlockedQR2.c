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
	int b= 2;
	printf(" b = %d \n", b);

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
for( int k =0 ; k< n; k++){

	/* Calculate number of blocks below diagonal */
	int nbloc = hn_bloc - k ; 
	int nbloc_copy = nbloc;

	/* Set initial c1 and c2, leftover constants */
	int c1 = 1; int c2 = 2; int leftover =0; 

	int iteration = 1;

	printf(" k = %d \n", k);		
	
	while( c2 <= nbloc ){
		
		for(int i=k; i< nbloc+k-1; i += c2){
		printf("k = %d, nbloc = %d, iteration = %d, Blocks called : %d, %d \n", k, nbloc, iteration, i, i+c1); 	
		}

		int mod = nbloc_copy % c2 ;
		nbloc_copy -= mod;		
		
		if ( mod != 0){
			if (leftover != 0){
			printf("k = %d, nbloc = %d, iteration = %d, Blocks called: %d, %d \n", k, nbloc, iteration, k+nbloc-c1, leftover);
			leftover =0; 	
			}else{
			leftover = k + nbloc -c1;
			}
		}

		/* Update constants c1 and c2 */
		c1 *= 2; 
		c2 *= 2;
		iteration += 1;

	}

	if (leftover != 0){
	printf("k = %d, nbloc = %d, Blocks called : %d, %d \n", k, nbloc, k, leftover); 
	leftover =0;
	}

	if( k == hn_bloc -1){
	printf(" k = %d, nbloc = %d, Blocks called : %d \n", k, nbloc, k);
	
	}	

		
}








}
