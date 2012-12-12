/* HPC 2012 Final Project: Jacqueline Bush, Paul Torres */

/* Functions that didn't really have a place anywhere else were placed here. */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "Utilities.h"


void prettyPrint(const double * const A, const int m, const int n){
  for(int i = 0; i < m ; i++){
    for(int j = 0; j < n; j++){
      printf("%5.5f ", A[i+j*m]);
    }
    printf("\n");
  }
  printf("\n");
}


/*---------------------Code to work with blocks of Matrix---------------------*/
 
void BlockMatrix(const double *inA, double *outA, int hA, int wA, int b, int i_bloc, int j_bloc){
/*----------------------------------------------------------------------------- 
PURPOSE: Takes matrix A, outputs the desired block, pads with zeros when necessary. 
ARGUEMENTS:
	hA: Height of matrix
	wA: Width of matrix  
	b: Blocksize
	i_bloc: i Block index 
	j_bloc: j Block index  
-----------------------------------------------------------------------------*/
int i, j;

	int wn_bloc = (wA+b-1)/b; // Number of Blocks 
	int hn_bloc = (hA+b -1)/b; // Number of Blocks
	int wpadding = wn_bloc*b - wA ; // Number of columns of zeros needed.
	int hpadding = hn_bloc*b - hA ; // Number of rows of zeros needed.
	for(i=0; i<b; i++){
		for(j=0; j<b; j++){
		
		if( (i_bloc == hn_bloc -1) && (i>= (b-hpadding))){
			outA[i+ j*b] = 0;
		} else{
			if( (j_bloc == wn_bloc -1) && (j >= (b-wpadding))){
			outA[i+ j*b] = 0;
			} else{
			outA[i + j*b] = inA[ i+i_bloc*b + (j+ j_bloc*b)*hA ] ;
			} 
		} }
    	} 


}


void UnBlockMatrix(double *outA, const double *inA, int hA, int wA, int b, int i_bloc, int j_bloc){
/*----------------------------------------------------------------------------- 
PURPOSE: Takes block of A, Replaces the values of A. Gets rid of pads if necessary.  
ARGUEMENTS:
	n: Size of matrix 
	b: Blocksize
	i_bloc: i Block index 
	j_bloc: j Block index  
-----------------------------------------------------------------------------*/

   int i, j;

        int wn_bloc = (wA+b -1)/b; // Number of Blocks 
	int hn_bloc = (hA+b -1)/b; // Number of Blocks
	int wpadding = wn_bloc*b - wA ; // Number of columns of zeros needed.
	int hpadding = hn_bloc*b - hA ; // Number of rows of zeros needed.

    if( ( wpadding == 0 && hpadding == 0) || (i_bloc != hn_bloc-1 && j_bloc != wn_bloc -1)){
	for(i=0; i<b; i++){
		for(j=0; j<b; j++){
		outA[ i+i_bloc*b + (j+j_bloc*b)*hA ] = inA[i + j*b];
		}
	} 
    }else{
	if( i_bloc == hn_bloc -1 && j_bloc != wn_bloc -1){
	for(i=0; i<(b-hpadding); i++){
		for(j=0; j<b; j++){
		outA[ i+i_bloc*b + (j+j_bloc*b)*hA ] = inA[i + j*b];
		}
	}
	}
	if( j_bloc == wn_bloc -1 && i_bloc != hn_bloc -1){
	for(i=0; i<b; i++){
		for(j=0; j<(b-wpadding); j++){
		outA[ i+i_bloc*b + (j+j_bloc*b)*hA ] = inA[i + j*b];
		}
	}
	}
	if( i_bloc == hn_bloc -1 && j_bloc == wn_bloc -1){
	for(i=0; i<(b-hpadding); i++){
		for(j=0; j<(b-wpadding); j++){
		outA[ i+i_bloc*b + (j+j_bloc*b)*hA ] = inA[i + j*b];
		}
	}
    	}
   }	
}

void CleanMatrix( double *A, int wA, int hA){
/*------------------------------------------------------------------------------
PURPOSE: Takes matrix A, Replaces it with zeros
ARGUEMENTS:
	wA: width of A
	hA: height of A
------------------------------------------------------------------------------*/

int i, j;

for(i=0; i<hA; i++){
	for(j=0; j<wA; j++){
	A[ i + j*hA] = 0 ;
	}
}

}


/*---------------------Code to work with blocks of Matrix---------------------*/
/*----------------------------------------------------------------------------- 
  PURPOSE: Takes vector A, outputs the desired block, pads with zeros when necessary. 
  ARGUMENTS:
  hA: Height of vector 
  b: Blocksize
  i_bloc: i Block index 
  -----------------------------------------------------------------------------*/
void BlockVectorForMatrixVector(const double *inA, double *outA, const int hA,
				const int b, const int i_bloc){
  int i;
  int hn_bloc = (hA+b-1)/b; // Number of Blocks
  int hpadding = hn_bloc*b - hA ; // Number of rows of zeros needed.
  //int hw_bloc = (wA+b-1)/b; // Number of Blocks
  //int wpadding = hn_bloc*b - wA ; // Number of rows of zeros needed.
  for(i = 0; i < b; i++){
    if((i_bloc == hn_bloc - 1) && (i >= (b - hpadding))){
      outA[i] = 0;
    } else {
      //do not bother padding for width...won't use those values anyway
      //omit all j-related terms
      outA[i] = inA[i + i_bloc*b];
    }
  }
}

/*----------------------------------------------------------------------------- 
PURPOSE: Takes block of vector A, replaces the values of A.
ARGUMENTS:
	n: Size of vector 
	b: Blocksize
	i_bloc: i Block index 
-----------------------------------------------------------------------------*/
void UnBlockVectorForMatrixVector(double *outA, const double *inA, const int hA,
				  const int b, const int i_bloc){
  int i;
  int hn_bloc = (hA + b - 1)/b; // Number of Blocks
  //int wpadding = b - 1 ; // Number of columns of zeros  to ignore
  int hpadding = hn_bloc*b - hA ; // Number of rows of zeros needed.
  if(hpadding > 0) printf("\n=unblock= padding: %d b: %d ibloc: %d\n", hpadding, b, i_bloc);
  if(hpadding == 0 || (i_bloc != hn_bloc-1)){
    for(i = 0; i < b; i++){
      outA[i + i_bloc*b] = inA[i];
    } 
  } else {
    if(i_bloc == hn_bloc - 1){
      for(i = 0; i < (b - hpadding); i++){
	outA[i + i_bloc*b] = inA[i];
      }
    }
  }
}
