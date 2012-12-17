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
#define verbose 1

void BlockedQR2( double *A, int h, int w, double *Q){

/* Block size */
	int b= 2;
	printf(" b = %d h = %d w = %d\n", b, h, w);
	
/* Calculate Number of blocks */
	int wn_bloc = (w+b-1)/b; // Number of Blocks in width (round up)
	int hn_bloc = (h+b -1)/b; // Number of Blocks in height (round up)
	//int * listofroots = malloc(sizeof(int)*20);
	//	for(int i = 0; i < 20; i++) 
	//listofroots[i] = 0;
/* Determine which num. of blocks is smaller */
	int n = wn_bloc; 
	if( wn_bloc > hn_bloc )
		n = hn_bloc;
	printf("n is %d \n", n);
/* Set Q = m by m identity matrix */
 	for(int i = 0; i<h; i++){
		for(int j=0; j<h; j++){
		Q[i + j*h] = 0;
		if(i == j)
			Q[i +j*h] = 1;
		}
	}	

/* Enter Loop */
for( int k =0 ; k < n; k++){
  /* Calculate number of blocks below diagonal */
  int nbloc = hn_bloc;
  /* Set initial c1 and c2 */
  int c = 1;
  int c1 = 1; 
  int c2 = 2; 
  int boundary = k;
  int iteration = 0;

  while(boundary < hn_bloc){
    c = 1;
    while(c < hn_bloc - boundary + 1){
      c*=2;
    }
    c/=2;
    c1 = 1;
    c2 = 2;
    printf(" x (nd/rd/st/th) boundary at %d, c: %d \n", boundary, c);
    printf(" k = %d, nbloc = %d \n", k, nbloc);		
    if(verbose) printf("c2 = %d and c: %d  hn_bloc: %d \n", c2, c, hn_bloc);
    while( c2 <= c ){
      printf("c2 = %d boundary = %d next boundary: %d\n", c2, boundary, boundary + c);
      for(int i = boundary; i < boundary + c - c1; i += c2){
	printf("k = %d, nbloc = %d, iteration = %d, BLOCKS called : %d, %d \n", k, nbloc, iteration, i, i+c1);
      }
      
      c1 *= 2; 
      c2 *= 2;
      iteration += 1;
    }
    if(boundary > 0){
      printf("\t k = %d, nbloc = %d, BLOCKS called : %d, %d \n", k, nbloc, k, boundary);    
    }
    boundary += c;
    if(boundary == hn_bloc - 1){
      printf("\tboundary = %d, nbloc = %d, BLOCKS called : %d \n", boundary, nbloc, boundary);    
    }
  }
  iteration = 0;

 }
}

