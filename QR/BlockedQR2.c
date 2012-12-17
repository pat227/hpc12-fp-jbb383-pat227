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
	int * listofroots = malloc(sizeof(int)*20);
	//	for(int i = 0; i < 20; i++) 
	//listofroots[i] = 0;
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
for( int k =0 ; k < 1; k++){
  /* Calculate number of blocks below diagonal */
  int nbloc = hn_bloc - k ; 
  /* Set initial c1 and c2 */
  int c = 2;
  int c1 = 1; 
  int c2 = 2; 
  int power = 0;
  int boundary = 0;
  int iteration = 0;
  int bigiter = 0;
  int localn = n - boundary - 1;
  int oldboundary = boundary;
  //handle special case of only 1 remaining block
  while(boundary <= n){
    if(localn == 1){
      printf("STOP HERE: k = %d, nbloc = %d, iteration = %d, Blocks called : %d, %d \n", k, nbloc, iteration, 0, oldboundary+1);
      break;
    }
    //compute how many we can evenly fit into binary tree and just do those
    while(c <= localn){
      c*=2;
      power++;
    }
    oldboundary = boundary;
    boundary += c;
    c1 = 1;
    c2 = 2;
    printf(" %d (rd/st/th) boundary at %d \n", bigiter, boundary);
    printf(" k = %d, nbloc = %d \n", k, nbloc);		
    
    while( c2 <= c ){
      printf("c2 = %d \n", c2);
      //only on first iteration do we start at zero, for rest it's prior boundary + 1
      if(bigiter > 0){
	for(int i=oldboundary+1; i < boundary; i += c2){
	  printf("k = %d, nbloc = %d, iteration = %d, Blocks called : %d, %d \n", k, nbloc, iteration, i, i+c1); 	
	}
      } else {
	for(int i=k; i < boundary; i += c2){
	  printf("k = %d, nbloc = %d, iteration = %d, Blocks called : %d, %d \n", k, nbloc, iteration, i, i+c1); 	
	}
      }
      
      /* Update constants c1 and c2 */
      c1 *= 2; 
      c2 *= 2;
      iteration += 1;
    }
    bigiter++;
    //now repeat for rest of "tree"
    iteration = 0;
    c = 2;
    localn = n - boundary;
    power = 0;
  }
 }
}
  /*
  while(c < localn){
    c*=2;
    power++;
  }
  int oldboundary = boundary;
  boundary += c;
  c2 = 2;
  c1 = 1;
  printf(" second boundary at %d \n", boundary);		
  printf(" k = %d, nbloc = %d \n", k, nbloc);		
	
  while( c2 <= c){
    printf("c2 = %d \n", c2);
    for(int i=oldboundary+1; i < boundary; i += c2){
      printf("k = %d, nbloc = %d, iteration = %d, Blocks called : %d, %d \n", k, nbloc, iteration, i, i+c1); 	
    }
    
    // Update constants c1 and c2
    c1 *= 2; 
    c2 *= 2;
    iteration += 1;
    
	}
  if(bigiter > 0){
    printf("k = %d, nbloc = %d, iteration = %d, Blocks called : %d, %d \n", k, nbloc, iteration, 0, oldboundary+1); 	
  }
  
 }








}
  */
