/* HPC 2012 Final Project: Jacqueline Bush, Paul Torres */

/* Functions that didn't really have a place anywhere else were placed here. */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "Utilities.h"
#include <math.h>

void prettyPrint(const double * const A, const int m, const int n){
  for(int i = 0; i < m ; i++){
    for(int j = 0; j < n; j++){
      printf("%5.5f ", A[i+j*m]);
    }
    printf("\n");
  }
  printf("\n");
}

//h is height NOT how many; how many is start row - h, computed by function
void copyColVector(double * source, double * target, int start, int h, int col){
  for(int i = start ; i<h; i++){
    target[i + col*h] = source[i];
  }
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
/*
  fname -> name of the file to which to write (append) the data
  m -> m size of matrix
  n -> n size of matrix
  iterations -> number of iterations for which algo was run (written in log-base-2 form)
  dependent -> dependent variable, could be time elapsed in seconds or gb/s
 */
void writetofile(const char * const fname, int m, int n, int iterations, double dependent){
  FILE * pf;
  char buffer[32];
  pf = fopen (fname,"a");
  double log = 0.0;
  if(pf!=NULL){
    //need a newline between series or else lines in gnu plot get screwed up
    fputs("\n", pf);
    //number of elements as log-base-2
    sprintf(buffer, "%f", (log10(m*n) / log10(2)) );
    fputs(buffer, pf);
    fputs(" ", pf);
    //use the log-base-10 of iterations
    log = log10(iterations);
    sprintf(buffer, "%f", log);
    fputs(buffer, pf);
    fputs(" ", pf);
    //the dependent variable - along z axis (typically time or gb/s)
    sprintf(buffer, "%f", dependent);
    fputs(buffer, pf);
    fputs("\n", pf);
    fclose(pf);
  } else {
    printf("Error opening file.");
    abort();
  }
}

void writetofile2(const char * const fname, int m, int n, double dependent){
  FILE * pf;
  char buffer[32];
  pf = fopen (fname,"a");
  double log = 0.0;
  if(pf!=NULL){
    //need a newline between series or else lines in gnu plot get screwed up
    fputs("\n", pf);
    sprintf(buffer, "%d", m);
    fputs(buffer, pf);
    fputs(" ", pf);
    sprintf(buffer, "%d", n);
    fputs(buffer, pf);
    fputs(" ", pf);
    //the dependent variable - along z axis (typically time or gb/s)
    sprintf(buffer, "%f", dependent);
    fputs(buffer, pf);
    fputs("\n", pf);
    fclose(pf);
  } else {
    printf("Error opening file.");
    abort();
  }
}

 /*----------------------------dgemm_simple Code-------------------------------*/
void dgemm_simple(const double *A, const int hA, const int wA, const double *B, const int hB, const int wB, double *C) {
/*----------------------------------------------------------------------------- 
PURPOSE: Computes simple matrix multiplication with A and B in Column-major order. 
ARGUEMENTS:
	wA: Width of A, number of columns in A
	hA: Height of A, number of rows in A
	wB: Width of B
	hB: Height of B 
-----------------------------------------------------------------------------*/
 
int hC = hA;
int wC = wB;

CleanMatrix(C , hC, wC);

for(int i=0; i<hC; i++){
	for(int j=0;j<wC; j++){
		for(int k=0; k<wA;k++){
			C[i + j*hC] += A[i +k*hA]*B[k +j*hB];
		}
	}
  } 
}

void simple_transpose(const double *A, int h, int w, double *B){
/*----------------------------------------------------------------------------- 
PURPOSE: Computes simple matrix transpose with A in Column-major order. 
ARGUEMENTS:
	w = width of A
	h = height of A
-----------------------------------------------------------------------------*/
  
  for(int j = 0; j < w; j++){
    for(int i = 0; i < h; i++){
      B[j + i*w] = A[i + j*h];
    }
  }
}
