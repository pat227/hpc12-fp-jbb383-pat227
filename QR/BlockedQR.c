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


void BlockedQR( double *A, int h, int w, double *Q, double *R){

int b = 3;

double *temp = malloc( sizeof(double)*2*b*2*b);
 

printf("A = \n");
prettyPrint(A, h, w);

CleanMatrix(temp, 2*b, 2*b);

UnBlockQ(temp, A, h, w, b, 0, 0, 1, 1 );



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

*/ 


void UnBlockQ( const double *inQ, double *outQ, int hQ, int wQ, int b, int i_bloc1, int j_bloc1, int i_bloc2, int j_bloc2){

printf("hi!\n");



}

void BlockQ(const double *inQ, double *outQ, int hQ, int wQ, int b, int i_bloc1, int j_bloc1, int i_bloc2, int j_bloc2 ){
/*----------------------------------------------------------------------------- 
PURPOSE: Takes matrix Q, outputs the desired 2b by 2b block, pads with zeros when necessary. 
ARGUEMENTS:
	hQ: Height of matrix
	wQ: Width of matrix  
	b: Blocksize
	The next four arguments tell us the two distinct columns and two distinct rows the new Q will effect. 
	i_bloc1: i Block index, uppermost i block index  
	j_bloc1: j Block index,  leftmost j block index
	i_bloc2: i Block index
	j_bloc2: j Block index 
-----------------------------------------------------------------------------*/

/* Check to make such i_bloc1 < i_bloc2 and j_bloc1 < j_block2 */
if (i_bloc1 >= i_bloc2)
  {
    fprintf(stderr, "We need i_bloc1 to be strictly greater than i_bloc2!\n");
    abort(); 
  } 


if (j_bloc1 >= j_bloc2)
 {
    fprintf(stderr, "We need j_bloc1 to be strictly greater than j_bloc2!\n");
    abort(); 
  } 


/* Calculate the number of rows or columns that need to be padded */
	int wn_bloc = (wQ+b-1)/b; // Number of Blocks 
	int hn_bloc = (hQ+b -1)/b; // Number of Blocks
	int wpadding = wn_bloc*b - wQ ; // Number of columns of zeros needed.
	int hpadding = hn_bloc*b - hQ ; // Number of rows of zeros needed.

	printf(" wn_bloc = %d \n", wn_bloc);
	printf(" hn_bloc = %d \n", hn_bloc);
	printf(" wpadding = %d \n", wpadding);
	printf(" hpadding = %d \n", hpadding);
	

/* Copy first Block - i_bloc1, j_bloc1 */
for(int i=0; i<b; i++){
		for(int j=0; j<b; j++){
		if( (i_bloc1 == hn_bloc -1) && (i>= (b-hpadding))){
			outQ[i + j*2*b] = 0;
		} else{
			if( (j_bloc2 == wn_bloc -1) && (j >= (b-wpadding))){
			outQ[i + j*2*b] = 0;
			} else{
			outQ[i + j*2*b] = inQ[ i+i_bloc1*b + (j+ j_bloc1*b)*hQ ] ;
			} 
		}	
		} 
	}

/* Copy second Block - i_bloc1, j_bloc2 */
for(int i=0; i<b; i++){
		for(int j=0; j<b; j++){
		if( (i_bloc1 == hn_bloc -1) && (i>= (b-hpadding))){
			outQ[2*b*b + i + j*2*b] = 0;
		} else{
			if( (j_bloc2 == wn_bloc -1) && (j >= (b-wpadding))){
			outQ[2*b*b + i + j*2*b] = 0;
			} else{
			outQ[2*b*b + i + j*2*b] = inQ[ i+i_bloc1*b + (j+ j_bloc2*b)*hQ ] ;
			} 
		}		
		} 
	}

/* Copy thrid Block - i_bloc2, j_bloc1 */
for(int i=0; i<b; i++){
		for(int j=0; j<b; j++){
		if( (i_bloc2 == hn_bloc -1) && (i>= (b-hpadding))){
			outQ[b+i + j*2*b] = 0;
		} else{
			if( (j_bloc1 == wn_bloc -1) && (j >= (b-wpadding))){
			outQ[b+i + j*2*b] = 0;
			} else{
			outQ[b+i + j*2*b] = inQ[ i+i_bloc2*b + (j+ j_bloc1*b)*hQ ] ;
			} 
		}
		} 
	}


/* Copy fourth Block - i_bloc2, j_bloc2 */
for(int i=0; i<b; i++){
		for(int j=0; j<b; j++){
		if( (i_bloc2 == hn_bloc -1) && (i>= (b-hpadding))){
			outQ[2*b*b +b +i + 2*j*b] = 0;
		} else{
			if( (j_bloc2 == wn_bloc -1) && (j >= (b-wpadding))){
			outQ[2*b*b +b +i + 2*j*b] = 0;
			} else{
			outQ[2*b*b +b +i + 2*j*b] = inQ[ i+i_bloc2*b + (j+ j_bloc2*b)*hQ ] ;
			} 
		}
		} 
	}
}






 









