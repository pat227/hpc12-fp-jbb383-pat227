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
  //int bigiter = 0;
  // int localhn = hn_bloc - boundary - k;
  //int oldboundary = boundary;
  //handle special case of only 1 remaining block
  while(boundary <= hn_bloc){
    //oldboundary = boundary;
    //if(verbose) printf("\n At top of loop: localhn is %d, oldboundary is %d", localhn, oldboundary);
    
    //if(localhn < 1){
      //printf("\nSTOP HERE A: k = %d, nbloc = %d, iteration = %d, Blocks called : %d, %d \n", k, nbloc, iteration, k, oldboundary);
    // break;
    //}
    
    //if(localhn < 2){
    //  printf("\nSTOP HERE B: k = %d, nbloc = %d, iteration = %d, Blocks called : %d, %d \n", k, nbloc, iteration, k, hn_bloc-1);
    //  break;
    //}
    //compute how many we can evenly fit into binary tree and just do those
    //might want to use math.c log function here & round down, if that would be faster
    c = 1;
    while(c <= hn_bloc - boundary + 1){
      c*=2;
    }
    
    c1 = 1;
    c2 = 2;
    printf(" x (nd/rd/st/th) boundary at %d \n", boundary);
    printf(" k = %d, nbloc = %d \n", k, nbloc);		
    if(verbose) printf("c2 = %d and c: %d  hn_bloc: %d \n", c2, c, hn_bloc);
    while( c2 <= c ){
      printf("c2 = %d \n", c2);
      //guard against going over the cliff
      //int loopguard = boundary;
      //if(boundary > hn_bloc) loopguard = hn_bloc;
      //      if(verbose) printf("loopguard is : %d\n", loopguard);
      //only on first iteration do we start at zero, for rest it's prior boundary + 1
      //if(bigiter > 0){
      //if(verbose) printf("Not first iteration...\n");
      for(int i = boundary; i < hn_bloc-c1; i += c2){
	printf("k = %d, nbloc = %d, iteration = %d, Blocks called : %d, %d \n", k, nbloc, iteration, i, i+c1);
      }
	//} else {
	//printf("First iteration...\n");
	//for(int i = k; i < loopguard; i += c2){
	//  printf("k = %d, nbloc = %d, iteration = %d, Blocks called : %d, %d \n", k, nbloc, iteration, i, i+c1);
	//	}
	//}
      
      /* Update constants c1 and c2 */
      c1 *= 2; 
      c2 *= 2;
      iteration += 1;
    }
    //call root of this binary tree against root of prior binary tree--which is always at zero
    //    if(bigiter>0){
    //printf("\nCall root of this most recent binary tree and root of all prior trees--always zero.");
    //  printf("\nk = %d, nbloc = %d, iteration = %d, Blocks called : %d, %d \n", k, nbloc, iteration, 0, oldboundary);
  boundary += c;
  }
  //bigiter++;
  //now repeat for rest of blocks that are own "trees"
  iteration = 0;
  //localhn = hn_bloc - boundary - k;
  //if(verbose) printf("\nBoundary is now (at end of while loop) % d and localhn is %d ", boundary);
  //printf("\nlocalhn is %d", localhn);
  if(k == hn_bloc - 1 ){
    printf("k = %d, nbloc = %d, Blocks called : %d \n", k, nbloc, k);    
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
