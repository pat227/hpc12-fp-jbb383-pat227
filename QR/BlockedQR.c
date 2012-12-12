
/* HPC 2012 Project  : Jacqueline Bush, Paul Torres */

/* BlockedQR FILE: 
Performs Blocked QR factorization/
*/

/* Headers: */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "WY.h"
#include "test.h"
#include "MatrixMatrixMultiply.h"
#include "MatrixTranspose.h"
#include "Utilities.h"
#include "BlockedQR.h"


void BlockedQR( double *A, int h, int w, double *Q){

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

/* Initialize Blocks */
	double *Q_block_1b_by_1b = malloc( b * b * sizeof( double));
	double *Q1_block_h_by_1b = malloc( h* b* sizeof(double));
	double *Q2_block_h_by_1b = malloc( h* b* sizeof(double));
	double *Q_block_2b_by_2b = malloc( 4*b * b * sizeof( double));
	double *Q1_block_h_by_2b = malloc( 2*h* b* sizeof(double));
	double *Q2_block_h_by_2b = malloc( 2*h* b* sizeof(double));
	double *Qt_block_1b_by_1b = malloc( b * b * sizeof( double));
	double *Qt_block_2b_by_2b = malloc( 2*b * 2 * b * sizeof(double));
	double *R_block_1b_by_1b = malloc( 2*b* b* sizeof(double));
	double *R_block_2b_by_1b = malloc(2 * b* b* sizeof(double)); 
	double *A_block_1b_by_1b = malloc(b * b* sizeof(double));
	double *A_block_2b_by_1b = malloc(2*b * b* sizeof(double));


/* Enter Loop */
for(int k =0; k< n; k++){	
	/* Update Diagonal Block */	
	BlockMatrix(A, A_block_1b_by_1b, h, w, b, k, k);
	WY( A_block_1b_by_1b, b, b, Q_block_1b_by_1b, Qt_block_1b_by_1b, R_block_1b_by_1b); 
	UnBlockMatrix(A, R_block_1b_by_1b, h, w, b, k,k);

	/* Update Corresponding Diagonal Block in in Q */
	BlockQ1(Q, Q1_block_h_by_1b, h, h, b, k );
	dgemm_simple(Q1_block_h_by_1b, h, b, Q_block_1b_by_1b, b, b, Q2_block_h_by_1b);
	//MatrixMatrixMultiply(Q1_block_h_by_1b, h, b, Q_block_1b_by_1b, b, b, Q2_block_h_by_1b);	
	UnBlockQ1(Q2_block_h_by_1b , Q, h, h, b, k);
			
	/* Update Blocks along row with Digonal Block */
	for( int j=k+1; j<wn_bloc ; j++){
	BlockMatrix(A, A_block_1b_by_1b, h, w, b, k, j);
	dgemm_simple(Qt_block_1b_by_1b, b, b, A_block_1b_by_1b, b, b, R_block_1b_by_1b);	
	//MatrixMatrixMultiply(Qt_block_1b_by_1b, b, b, A_block_1b_by_1b, b, b, R_block_1b_by_1b);
	UnBlockMatrix(A, R_block_1b_by_1b, h,w,b, k, j);
	}

	/* Update Blocks below the Diagonal (will also update diagonal) */
	for( int i=k+1; i<hn_bloc ; i++){
		Block(A_block_2b_by_1b, A, h , w, b, k, i, k);	
		WY( A_block_2b_by_1b, 2*b, b, Q_block_2b_by_2b, Qt_block_2b_by_2b, R_block_2b_by_1b);
		UnBlock(A, R_block_2b_by_1b, h, w, b,k, i, k);

		/* Update Q */
		BlockQ(Q, Q1_block_h_by_2b, h, h, b, k, i );
		dgemm_simple(Q1_block_h_by_1b, 2*b, h, Q_block_2b_by_2b, 2*b, 2*b, Q2_block_h_by_2b);
		//MatrixMatrixMultiply(Q1_block_h_by_2b, h, 2*b, Q_block_2b_by_2b, 2*b, 2*b, Q2_block_h_by_2b);		
		UnBlockQ( Q2_block_h_by_2b, Q, h, h, b, k, i);	
 
		/* Update Blocks along i and k  rows */
		for( int j=k+1; j<wn_bloc ; j++){
			Block(A_block_2b_by_1b, A , h, w, b, k, i, j);
			dgemm_simple(Qt_block_2b_by_2b, 2*b, 2*b, A_block_2b_by_1b, 2*b, b, R_block_2b_by_1b);
			//MatrixMatrixMultiply(Qt_block_2b_by_2b, 2*b, 2*b, A_block_2b_by_1b, 2*b, b, R_block_2b_by_1b);
			UnBlock(A, R_block_2b_by_1b, h, w, b,k, i, j);	
		}

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






 









