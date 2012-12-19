
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


void BlockedQR( double *A, int h, int w, double *Q){

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


/* Enter Loop */
for(int k =0; k< n; k++){	
		
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
	if( (k == hn_bloc - 1) && (n == hn_bloc) ){
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
	}
}

}

/*---------------------Code to work with blocks of Matrix---------------------*/

/* Note: Because more than one code in this project blocks and unblocks an individual block of a matrix, the code to block and unblock ONE block of the matrix is in the utilities file.  In this section we deal with blocking issues particular to this code. 

The code in the utilities Folder is:
CleanMatrix( double *A, int wA, int hA)
--- Zeros out Matrix A

UnBlockMatrix(double *outA, const double *inA, int hA, int wA, int b, int i_bloc, int j_bloc)
--- Copies elements in block into appropriate spot in matrix

BlockMatrix(const double *inA, double *outA, int hA, int wA, int b, int i_bloc, int j_bloc)
--- Copies appropriate spot in matrix into a block 

Here we have:
Block( double *outA, const double *inA, int h, int w, int b, int i_bloc1, int i_bloc2, int j_bloc)
--- Copies appropiate 2 spots in matrix into one block

UnBlock( double *outA, const double *inA, int h, int w, int b, int i_bloc1, int i_bloc2, int j_bloc)
--- Copies elements in block into appropriate 2 spots in matrix

BlockQ(const double *inQ, double *outQ, int hQ, int wQ, int b, int j_bloc1, int j_bloc2 )
--- Copies appropriate 2 sets of columns in matrix into one block

UnBlockQ( const double *inQ, double *outQ, int hQ, int wQ, int b, int j_bloc1, int j_bloc2)
--- Copies elements in block into appropriate spots in matrix 

BlockQ1(const double *inQ, double *outQ, int hQ, int wQ, int b, int j_bloc )
--- Copies appropriate set of columns in matrix into one block

UnBlockQ1( const double *inQ, double *outQ, int hQ, int wQ, int b, int j_bloc)
--- Copies elements in block into appropriate spots in matrix 


*/ 

void UnBlock( double *outA, const double *inA, int h, int w, int b, int i_bloc1, int i_bloc2, int j_bloc){
/*----------------------------------------------------------------------------- 
PURPOSE: Takes 2b by b block, Copies block into appropriate place in matrix, removes padding as necessary.  
ARGUEMENTS:
	h: Height of matrix
	w: Width of matrix  
	b: Blocksize
	The next four arguments tell us the two distinct columns and two distinct rows the new Q will effect. 
	i_bloc1: i Block index, uppermost i block index  
	i_bloc2: i Block index
	j_bloc: j Block index i.e. what column we are working in.
-----------------------------------------------------------------------------*/


/* Check to make such i_bloc1 < i_bloc2 and j_bloc1 < j_block2 */
if (i_bloc1 >= i_bloc2)
  {
    fprintf(stderr, "We need i_bloc1 to be strictly greater than i_bloc2!\n");
    abort(); 
  } 

/* Calculate the number of rows or columns padded that needs to be removed*/
	int wn_bloc = (w+b-1)/b; // Number of Blocks 
	int hn_bloc = (h+b -1)/b; // Number of Blocks
	int wpadding = wn_bloc*b - w ; // Number of columns of zeros needed.
	int hpadding = hn_bloc*b - h ; // Number of rows of zeros needed.


/* Check to make sure that i_bloc2 < hn_bloc and j_bloc < wn_bloc */
if (i_bloc2 >= hn_bloc)
  {
    fprintf(stderr, "We need i_bloc2 to be strictly less than the number of i blocks, %d!\n", hn_bloc);
    abort(); 
  } 

if (j_bloc >= wn_bloc)
 {
    fprintf(stderr, "We need j_bloc to be strictly less than the number of j blocks, %d!\n", wn_bloc);
    abort(); 
  } 


/* First Block - i_bloc1, j_bloc */
    if( ( wpadding == 0 && hpadding == 0) || (i_bloc1 != hn_bloc-1 && j_bloc != wn_bloc -1)){
	for(int i=0; i<b; i++){
		for(int j=0; j<b; j++){
		outA[ i+i_bloc1*b + (j+j_bloc*b)*h ] = inA[i + 2*j*b];
		}
	} 
    }else{
	if( i_bloc1 == hn_bloc -1 && j_bloc != wn_bloc -1){
	for(int i=0; i<(b-hpadding); i++){
		for(int j=0; j<b; j++){
		outA[ i+i_bloc1*b + (j+j_bloc*b)*h ] = inA[i + 2*j*b];
		}
	}
	}
	if( j_bloc == wn_bloc -1 && i_bloc1 != hn_bloc -1){
	for(int i=0; i<b; i++){
		for(int j=0; j<(b-wpadding); j++){
		outA[ i+i_bloc1*b + (j+j_bloc*b)*h ] = inA[i + 2*j*b];
		}
	}
	}
	if( i_bloc1 == hn_bloc -1 && j_bloc == wn_bloc -1){
	for(int i=0; i<(b-hpadding); i++){
		for(int j=0; j<(b-wpadding); j++){
		outA[ i+i_bloc1*b + (j+j_bloc*b)*h ] = inA[i + 2*j*b];
		}
	}
    	}
   }


/* Second Block, i_bloc2, j_bloc */
if( ( wpadding == 0 && hpadding == 0) || (i_bloc2 != hn_bloc-1 && j_bloc != wn_bloc -1)){
	for(int i=0; i<b; i++){
		for(int j=0; j<b; j++){
		outA[ i+i_bloc2*b + (j+j_bloc*b)*h ] = inA[b+i + j*2*b];
		}
	} 
    }else{
	if( i_bloc2 == hn_bloc -1 && j_bloc != wn_bloc -1){
	for(int i=0; i<(b-hpadding); i++){
		for(int j=0; j<b; j++){
		outA[ i+i_bloc2*b + (j+j_bloc*b)*h ] = inA[b+i + j*2*b];
		}
	}
	}
	if( j_bloc == wn_bloc -1 && i_bloc2 != hn_bloc -1){
	for(int i=0; i<b; i++){
		for(int j=0; j<(b-wpadding); j++){
		outA[ i+i_bloc2*b + (j+j_bloc*b)*h ] = inA[b+i + j*2*b];
		}
	}
	}
	if( i_bloc2 == hn_bloc -1 && j_bloc == wn_bloc -1){
	for(int i=0; i<(b-hpadding); i++){
		for(int j=0; j<(b-wpadding); j++){
		outA[ i+i_bloc2*b + (j+j_bloc*b)*h ] = inA[b+i + j*2*b];
		}
	}
    	}
   }


}


void Block( double *outA, const double *inA, int h, int w, int b, int i_bloc1, int i_bloc2, int j_bloc){
/*----------------------------------------------------------------------------- 
PURPOSE: Takes matrix A, outputs the desired 2b by b block, pads with zeros when necessary. 
ARGUEMENTS:
	h: Height of matrix
	w: Width of matrix  
	b: Blocksize
	The next four arguments tell us the two distinct columns and two distinct rows the new Q will effect. 
	i_bloc1: i Block index, uppermost i block index  
	i_bloc2: i Block index
	j_bloc: j Block index i.e. what column we are working in.
-----------------------------------------------------------------------------*/

/* Check to make such i_bloc1 < i_bloc2 and j_bloc1 < j_block2 */
if (i_bloc1 >= i_bloc2)
  {
    fprintf(stderr, "We need i_bloc1 to be strictly greater than i_bloc2!\n");
    abort(); 
  } 


/* Calculate the number of rows or columns that need to be padded */
	int wn_bloc = (w+b-1)/b; // Number of Blocks 
	int hn_bloc = (h+b -1)/b; // Number of Blocks
	int wpadding = wn_bloc*b - w ; // Number of columns of zeros needed.
	int hpadding = hn_bloc*b - h ; // Number of rows of zeros needed.

/* Copy first Block - i_bloc1, j_bloc */
for(int i=0; i<b; i++){
		for(int j=0; j<b; j++){
		if( (i_bloc1 == hn_bloc -1) && (i>= (b-hpadding))){
			outA[i + j*2*b] = 0;
		} else{
			if( (j_bloc == wn_bloc -1) && (j >= (b-wpadding))){
			outA[i + j*2*b] = 0;
			} else{
			outA[i + j*2*b] = inA[ i+i_bloc1*b + (j+ j_bloc*b)*h ] ;
			} 
		}	
		} 
	}


/* Copy second Block - i_bloc2, j_bloc */
for(int i=0; i<b; i++){
		for(int j=0; j<b; j++){
		if( (i_bloc2 == hn_bloc -1) && (i>= (b-hpadding))){
			outA[b+i + j*2*b] = 0;
		} else{
			if( (j_bloc == wn_bloc -1) && (j >= (b-wpadding))){
			outA[b+i + j*2*b] = 0;
			} else{
			outA[b+i + j*2*b] = inA[ i+i_bloc2*b + (j+ j_bloc*b)*h ] ;
			} 
		}
		} 
	}


}


void BlockQ1(const double *inQ, double *outQ, int hQ, int wQ, int b, int j_bloc ){
/*----------------------------------------------------------------------------- 
PURPOSE: Takes matrix Q, outputs the desired h by b block pads with zeros when necessary. 
ARGUEMENTS:
	hQ: Height of matrix
	wQ: Width of matrix  
	b: Blocksize 
	j_bloc: j Block index
-----------------------------------------------------------------------------*/


/* Calculate the number of rows or columns that need to be padded */
	int wn_bloc = (wQ+b-1)/b; // Number of Blocks 
	int wpadding = wn_bloc*b - wQ ; // Number of columns of zeros needed.



if (j_bloc >= wn_bloc)
 {
    fprintf(stderr, "BlockQ1: We need j_bloc2 to be strictly less than the number of j blocks, %d!\n", wn_bloc);
    abort(); 
  } 

/* Copy Block - j_bloc */
        for(int i=0; i<hQ; i++){
		for(int j=0; j<b; j++){
			if( (j_bloc == wn_bloc -1) && (j >= (b-wpadding))){
			outQ[i + j*hQ] = 0;
			} else{
			outQ[i + j*hQ] = inQ[ i + (j+ j_bloc*b)*hQ ] ;
			} 
		}	 
	}


}


void UnBlockQ1( const double *inQ, double  *outQ, int hQ, int wQ, int b, int j_bloc){
/*----------------------------------------------------------------------------- 
PURPOSE: Takes block of Q, Replaces the values of A. Gets rid of pads if necessary.  
ARGUEMENTS: 
	b: Blocksize
	hQ: height of Q
	wQ: width of Q  
	j_bloc: j Block index
-----------------------------------------------------------------------------*/


/* Calculate the number of rows or columns padded that needs to be removed*/
	int wn_bloc = (wQ+b-1)/b; // Number of Blocks 
	int wpadding = wn_bloc*b - wQ ; // Number of columns of zeros needed.


/* Check to make sure that j_bloc < wn_bloc */

if (j_bloc >= wn_bloc)
 {
    fprintf(stderr, "We need j_bloc2 to be strictly less than the number of j blocks, %d!\n", wn_bloc);
    abort(); 
  } 

/* Block - j_bloc */

    if( ( wpadding == 0) || ( j_bloc != wn_bloc -1)){
	for(int i=0; i<hQ; i++){
		for(int j=0; j<b; j++){
		outQ[ i + (j+j_bloc*b)*hQ ] = inQ[i + j*hQ];
		}
	} 
    }else{
	if( j_bloc == wn_bloc -1){
	for(int i=0; i<hQ; i++){
		for(int j=0; j<(b-wpadding); j++){
		outQ[ i+ + (j+j_bloc*b)*hQ ] = inQ[ i + j*hQ];
		}
	}
	}
	
     }



}

void UnBlockQ( const double *inQ, double *outQ, int hQ, int wQ, int b, int j_bloc1, int j_bloc2){
/*----------------------------------------------------------------------------- 
PURPOSE: Takes block of Q, Replaces the values of A. Gets rid of pads if necessary.  
ARGUEMENTS:
	n: Size of matrix 
	b: Blocksize
	hQ: height of Q
	wQ: width of Q  
	j_bloc1: j Block index,  leftmost j block index
	j_bloc2: j Block index 
-----------------------------------------------------------------------------*/

/* Check to make such j_bloc1 < j_block2 */

if (j_bloc1 >= j_bloc2)
 {
    fprintf(stderr, "We need j_bloc1 to be strictly greater than j_bloc2!\n");
    abort(); 
  } 


/* Calculate the number of rows or columns padded that needs to be removed*/
	int wn_bloc = (wQ+b-1)/b; // Number of Blocks 
	int wpadding = wn_bloc*b - wQ ; // Number of columns of zeros needed.


/* Check to make sure that j_bloc < wn_bloc */

if (j_bloc2 >= wn_bloc)
 {
    fprintf(stderr, "We need j_bloc2 to be strictly less than the number of j blocks, %d!\n", wn_bloc);
    abort(); 
  } 


/* First Block - j_bloc1 */
	for(int i=0; i<hQ; i++){
		for(int j=0; j<b; j++){
		outQ[i + (j+ j_bloc1*b)*hQ] = inQ[i + j*hQ];
		}
	}

/* Second Block - j_bloc2 */

    if( ( wpadding == 0) || ( j_bloc2 != wn_bloc -1)){
	for(int i=0; i<hQ; i++){
		for(int j=0; j<b; j++){
		outQ[ i + (j+j_bloc2*b)*hQ ] = inQ[hQ*b + i + j*hQ];
		}
	} 
    }else{
	if( j_bloc2 == wn_bloc -1){
	for(int i=0; i<hQ; i++){
		for(int j=0; j<(b-wpadding); j++){
		outQ[ i+ + (j+j_bloc2*b)*hQ ] = inQ[hQ*b + i + j*hQ];
		}
	}
	}
	
     }



}

void BlockQ(const double *inQ, double *outQ, int hQ, int wQ, int b, int j_bloc1, int j_bloc2 ){
/*----------------------------------------------------------------------------- 
PURPOSE: Takes matrix Q, outputs the desired 2b by 2b block, pads with zeros when necessary. 
ARGUEMENTS:
	hQ: Height of matrix
	wQ: Width of matrix  
	b: Blocksize 
	j_bloc1: j Block index,  leftmost j block index
	j_bloc2: j Block index 
-----------------------------------------------------------------------------*/

/* Check to make such j_bloc1 < j_block2 < wn_bloc.*/

if (j_bloc1 >= j_bloc2)
 {
    fprintf(stderr, "We need j_bloc1 to be strictly greater than j_bloc2!\n");
    abort(); 
  } 


/* Calculate the number of rows or columns that need to be padded */
	int wn_bloc = (wQ+b-1)/b; // Number of Blocks 
	int wpadding = wn_bloc*b - wQ ; // Number of columns of zeros needed.

if (j_bloc2 >= wn_bloc)
 {
    fprintf(stderr, "We need j_bloc2 to be strictly less than the number of j blocks, %d!\n", wn_bloc);
    abort(); 
  } 

	
/* Copy first Block - j_bloc1 */
      for(int i=0; i<hQ; i++){
		for(int j=0; j<b; j++){
			outQ[i + j*hQ] = inQ[ i + (j+ j_bloc1*b)*hQ ] ;
		}	 
	}

/* Copy second Block - j_bloc2 */
        for(int i=0; i<hQ; i++){
		for(int j=0; j<b; j++){
			if( (j_bloc2 == wn_bloc -1) && (j >= (b-wpadding))){
			outQ[hQ*b + i + j*hQ] = 0;
			} else{
			outQ[hQ*b + i + j*hQ] = inQ[ i + (j+ j_bloc2*b)*hQ ] ;
			} 
		}	 
	}


}






 









