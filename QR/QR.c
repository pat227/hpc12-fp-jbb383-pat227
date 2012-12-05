#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "QR.h"

#define L1_BLK_SIZE 2
#define MULTIPLE 2
#define L2_BLK_SIZE (L1_BLK_SIZE * MULTIPLE)

/* HPC 2012 Project : Jacqueline Bush, Paul Torres */
/*==================================================================================*/
/* Code Preforms Blocked Matrix Matrix Multiplication, Matrices do not need to be square */
/*==================================================================================*/
/* Blocked Matrix Matrix Multiply Subfunctions */
static void dgemm_lowest( const double*restrict, const double*restrict, double*restrict);
static void dgemm_middle( const double*restrict, const double*restrict, double*restrict); 

/* Code to load and unload blocks of matrix */
void BlockMatrix( const double *, double *, int, int, int, int, int);
void BlockMatrixInner(const double *inA, double *outA, int hA, int wA, int b, int i_bloc, int j_bloc);
void UnBlockMatrix( double *, const double *, int, int, int, int, int);
void CleanMatrix( double *, int, int );

/* Code to test Blocked Matrix Matrix Multipy */
void dgemm_simple(const double *, const int, const int, const double *, const int, const int, double *);
void test( const double *, const double*, double*, int, int, int , int );

int MatrixMatrixMultiply( double *A, int wA, int hA, double *B, int wB, int hB, double *C)
{

/* Check that Matrix deminsions are valid */
if( wA !=  hB){
{
    fprintf(stderr, "Matrix Multiplication Invalid! \nThe number of columns of A must equal the number of rows of B!\n");
    abort(); 
  }
}

/* Matrix Deminsions of Output matrix */
int wC = wB;
int hC = hA;

/* Calculate the number of Blocks in the height and width of C, as well as the width of A */
int wn_bloc = (wC+L2_BLK_SIZE -1)/L2_BLK_SIZE; // Number of Blocks in the widith (round up)
int hn_bloc = (hC+L2_BLK_SIZE -1)/L2_BLK_SIZE; // Number of Blocks in the height (round up)
int wA_bloc = (wA+L2_BLK_SIZE -1)/L2_BLK_SIZE; // Number of Blocks in the widith of matrix A  (round up)

int b = L2_BLK_SIZE; // Block Size

/* Make static matrices */
static __attribute__ ((aligned(16))) double a_block[L2_BLK_SIZE*L2_BLK_SIZE];
static __attribute__ ((aligned(16))) double b_block[L2_BLK_SIZE*L2_BLK_SIZE];
static __attribute__ ((aligned(16))) double c_block[L2_BLK_SIZE*L2_BLK_SIZE]; 

/* Clean matrices */
CleanMatrix(c_block, b, b);
CleanMatrix(a_block, b, b);
CleanMatrix(b_block, b, b);

for( int i = 0; i < hn_bloc; i++){
	for( int j=0; j < wn_bloc; j++){
		BlockMatrix(C, c_block, hC, wC,  b, i, j);
		for( int k = 0; k < wA_bloc ; k++){
		BlockMatrix(A, a_block, hA, wA, b, i, k);
		BlockMatrix(B, b_block, hB, wB, b, k, j);
		dgemm_middle(a_block, b_block, c_block);
		CleanMatrix(a_block, b, b);
		CleanMatrix(b_block, b, b);	
		}
		UnBlockMatrix(C, c_block, hC, wC, b, i,j);
		CleanMatrix(c_block, b, b);		
	}
}


/* Uncomment to Test this function */
//test(A, B, C, hA, wA, hB, wB);

return 0;
}


/*=============================================================================*/
 

/*----------------------------dgemm_middle Code-------------------------------*/
static void dgemm_middle(const double*restrict A, const double*restrict B, double*restrict C){
/*----------------------------------------------------------------------------- 
PURPOSE: Takes fixed sized matrices A and B, size L2_BLK_SIZE, Outputs result matrix C.  
ARGUEMENTS:
-----------------------------------------------------------------------------*/

int n = L2_BLK_SIZE; // Size of Matrix 
int n_bloc = L2_BLK_SIZE/L1_BLK_SIZE; // Number of Blocks
int b = L1_BLK_SIZE; // Block Size

/* Make static matrices */
static __attribute__ ((aligned(16))) double a_block[L1_BLK_SIZE*L1_BLK_SIZE];
static __attribute__ ((aligned(16))) double b_block[L1_BLK_SIZE*L1_BLK_SIZE];
static __attribute__ ((aligned(16))) double c_block[L1_BLK_SIZE*L1_BLK_SIZE]; 

CleanMatrix(c_block, b, b);

for( int j = 0; j < n_bloc; j++){
	for( int i=0; i < n_bloc; i++){
		BlockMatrixInner(C, c_block, n, n, b, i, j);
		for( int k = 0; k < n_bloc ; k++){
		BlockMatrixInner(A, a_block, n, n, b, i, k);
		BlockMatrixInner(B, b_block, n, n, b, k, j);
		dgemm_lowest( a_block, b_block, c_block);
		CleanMatrix(a_block, b, b);
		CleanMatrix(b_block, b, b);	
		}
		UnBlockMatrix(C, c_block, n, n, b, i,j);		
	}
}	

}


/*----------------------------dgemm_lowest Code-------------------------------*/
static void dgemm_lowest(const double*restrict A, const double*restrict B, double*restrict C){
/*----------------------------------------------------------------------------- 
PURPOSE: Takes fixed sized matrices A and B, size L1_BLK_SIZE, Outputs result matrix C.  
ARGUEMENTS:
-----------------------------------------------------------------------------*/

int n = L1_BLK_SIZE , i, j, k;
 
  for(j=0; j<n; j++){
	for(k=0;k<n; k++){
		for(i=0; i<n;i++){
		C[i + j*n] += A[i + k*n]*B[k + j*n];
		}
	}
  }

}

/*============================================================================*/


/*----------------------------dgemm_simple Code-------------------------------*/
void dgemm_simple(const double *A, const int wA, const int hA, const double *B, const int wB, const int hB, double *C) {
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

for(int i=0; i<hC; i++){
	for(int j=0;j<wC; j++){
		for(int k=0; k<wA;k++){
			C[i + j*hC] += A[i +k*hA]*B[k +j*hB];
		}
	}
  }
  
}

void test( const double *A, const double *B, double *C, int hA, int wA, int hB, int wB ){
/*----------------------------------------------------------------------------- 
PURPOSE: Tests Blocked matrix multiplication against the simple variant (which we know works)
ARGUEMENTS:
	hA: Height of A
	wA: Width of A
	hB: Height of B
	wB: Width of B
-----------------------------------------------------------------------------*/

/* Dimensions of C */
int hC = hA;
int wC = wB;

/* Test Matrix */
double *testC = malloc( wC*hC*sizeof(double) );
CleanMatrix(testC, hC, wC);
dgemm_simple( A, wA, hA, B, wB, hB, testC);

		
/* Error Check */
for(int i=0; i<hC; i++){
	for(int j=0;j<wC; j++){
		double error = abs( C[i+ j*hC] - testC[i + j*hC]);
		//printf(" error = %f, i = %d, j = %d\n", error, i, j);
		//printf(" C = %f, testC = %f \n", C[i+j*hC], testC[i + j*hC]);
		double errorbound = 1e-5;
		if( error > errorbound ){
		fprintf(stderr,"Blocked Matrix Multiplication is not working! \n") ;
   	 	abort();
		}
	}
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

if( wA % b != 0 || hA % b != 0  ){
	// Will only enter this section if we are in outerblock of code. 
	int wn_bloc = (wA+L2_BLK_SIZE -1)/L2_BLK_SIZE; // Number of Blocks 
	int hn_bloc = (hA+L2_BLK_SIZE -1)/L2_BLK_SIZE; // Number of Blocks
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
}else{
	for(i=0; i<b; i++){
		for(j=0; j<b; j++){
		outA[i + j*b] = inA[ i+i_bloc*b + (j+ j_bloc*b)*hA ] ;
		} 
	}
}

}

void BlockMatrixInner(const double *inA, double *outA, int hA, int wA, int b, int i_bloc, int j_bloc){
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

if( wA % b != 0 || hA % b != 0  ){
	// Will only enter this section if we are in outerblock of code. 
	int wn_bloc = (wA+L2_BLK_SIZE -1)/L2_BLK_SIZE; // Number of Blocks 
	int hn_bloc = (hA+L2_BLK_SIZE -1)/L2_BLK_SIZE; // Number of Blocks
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
}else{
	for(i=0; i<b; i++){
		for(j=0; j<b; j++){
		outA[i + j*b] = inA[ i+i_bloc*b + (j+ j_bloc*b)*hA ] ;
		} 
	}
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



//==============================================================================
// Matrix transpose and associated functions
//==============================================================================
void copyTransposedL2Block(const double * const B, const int submatrix_row , 
			   const int submatrix_col, double * A, const int w, 
			   const int h);
void copyL2Block(const double * const A, const int submatrix_row , 
		 const int submatrix_col, double * B, const int w, 
		 const int h);
void copyL1Block(const double * const A, const int i, const int j, double * B);
void transposeL2(const double * const A, double * B);
void copyTransposedL1Block(const double * const B, const int submatrix_row, 
			   const int submatrix_col, double * A);
void transposeL1Size(double * const A, double * B);
void prettyPrint(const double * const A, const int m, const int n){
  //printf("m:%d n:%d \n",m,n);
  for(int i = 0; i < m ; i++){
    for(int j = 0; j < n; j++){
      printf("%5.5f ", A[i+j*m]);
    }
    printf("\n");
  }
  printf("\n");
}
/*
  Transposes a matrix using L1 and L2 sized sub-blocks to hopefully gain a speed
  up, similar to the openCL matrix transposition although that was on a GPU...
  Matrices can be square or not.
*/

//Answer matrix B must have the dimensions of A transposed--function does not 
//check for it; function transposes the matrix by blocks; each element of each
//submatrix block is transposed and then each block is copied back to its 
//transposed location within the larger matrix
void MatrixTranspose(const double * const A, const int h, const int w, double * B){
  int num_blocks_w = (w + L2_BLK_SIZE - 1) / L2_BLK_SIZE;
  int num_blocks_h = (h + L2_BLK_SIZE - 1) / L2_BLK_SIZE;
  int verbose = 0;
  if(verbose) printf("\n# of sub-blocks:  %d high by %d wide", num_blocks_h, num_blocks_w);
  /* Make static matrices */
  static __attribute__ ((aligned(16))) double b_blockL2[L2_BLK_SIZE*L2_BLK_SIZE];
  static __attribute__ ((aligned(16))) double c_blockL2[L2_BLK_SIZE*L2_BLK_SIZE];
  for(int j = 0; j < num_blocks_w; j++){
    for(int i = 0; i < num_blocks_h; i++){
      copyL2Block(A, i , j, b_blockL2, w, h);
      if(verbose) {
	printf("\nb_blockL2:\n");
	prettyPrint(b_blockL2, L2_BLK_SIZE, L2_BLK_SIZE);
	printf("\nInvoking tranposeL2");
      }
      transposeL2(b_blockL2, c_blockL2);
      if(verbose){
	printf("\nc_blockL2:\n");
	prettyPrint(c_blockL2, L2_BLK_SIZE, L2_BLK_SIZE);
      }
      copyTransposedL2Block(c_blockL2, i, j, B, w, h);
      if(verbose){
	printf("\nUpdated B:\n");
	prettyPrint(B, w, h);
      }
    }
  }
}

//copy a transposed L2-sized block to proper location in larger answer matrix--
//submatrices do not go back to where they came from but are themselves 
//transposed as well
void copyTransposedL2Block(const double * const B, const int submatrix_row , 
			   const int submatrix_col, double * A, const int w, 
			   const int h){
  int verbose = 0;
  if(verbose){
    printf("\n=copyTransposedL2Block=");
    printf("\nsubmatrix_col: %d, h: %d L2: %d, w: %d submatrix_row: %d", 
	   submatrix_col, h, L2_BLK_SIZE, w, submatrix_row);
  }
  if((submatrix_col+1)*L2_BLK_SIZE <= w && (submatrix_row+1)*L2_BLK_SIZE <= h){
    if(verbose) printf("=copyTransposedL2Block= case 1");
    for(int j = 0; j < L2_BLK_SIZE; j++){
      for(int i = 0; i < L2_BLK_SIZE; i++){
	A[i + (submatrix_col*L2_BLK_SIZE) + (submatrix_row*L2_BLK_SIZE+j)*w] = 
	  B[i+j*L2_BLK_SIZE];
      }
    } 
  } else if((submatrix_col+1)*L2_BLK_SIZE > w && 
	    (submatrix_row+1)*L2_BLK_SIZE <= h){
    if(verbose) printf("\n=copyTransposedL2Block= case 2");
    //have to do stop early along bottom edge at some point
    int diff = (submatrix_col+1)*L2_BLK_SIZE - w;
    int new_rowguard = L2_BLK_SIZE - diff;
    for(int j = 0; j < L2_BLK_SIZE; j++){
      for(int i = 0; i < new_rowguard; i++){
	A[i + (submatrix_col*L2_BLK_SIZE) + (submatrix_row*L2_BLK_SIZE+j)*w] = 
	  B[i+j*L2_BLK_SIZE];
      }
    } 
  } else if((submatrix_col+1)*L2_BLK_SIZE <= w && 
	    (submatrix_row+1)*L2_BLK_SIZE > h){
    if(verbose) printf("\n=copyTransposedL2Block= case 3");
    int diff = (submatrix_row+1)*L2_BLK_SIZE - h;
    int new_colguard = h - diff;
    for(int j = 0; j < new_colguard; j++){
      for(int i = 0; i < L2_BLK_SIZE; i++){
	A[i + (submatrix_col)*L2_BLK_SIZE + (submatrix_row*L2_BLK_SIZE+j)*w] = 
	  B[i+j*L2_BLK_SIZE];
      }
    } 
  } else if((submatrix_col+1)*L2_BLK_SIZE > w && 
	    (submatrix_row+1)*L2_BLK_SIZE > h){
    if(verbose) printf("\n=copyTransposedL2Block= case 4");
    int diff1 = (submatrix_col+1)*L2_BLK_SIZE - h;
    int new_rowguard = L2_BLK_SIZE - diff1;
    int diff2 = (submatrix_row+1)*L2_BLK_SIZE - w;
    int new_colguard = w - diff2;
    for(int j = 0; j < new_colguard; j++){
      for(int i = 0; i < new_rowguard; i++){
	A[i + (submatrix_col*L2_BLK_SIZE) + (submatrix_row*L2_BLK_SIZE+j)*w] =
	  B[i+j*L2_BLK_SIZE];
      }
    } 
  }
}
//copy a NOT-YET-transposed L2 submatrix block of A to L2-block-sized B
void copyL2Block(const double * const A, const int submatrix_row , 
		 const int submatrix_col, double * B, const int w, const int h){
  int verbose = 0;
  if(verbose){
    printf("\n=copyL2Block=");
    printf("\nsubmatrix_col: %d, h: %d L2: %d, w: %d submatrix_row: %d", 
	   submatrix_col, h, L2_BLK_SIZE, w, submatrix_row);
  }
  if((submatrix_row+1)*L2_BLK_SIZE <= h && (submatrix_col+1)*L2_BLK_SIZE <= w){
    if(verbose) printf("\n=copyL2Block= case 1");
    for(int j = 0; j < L2_BLK_SIZE ; j++){
      for(int i = 0; i < L2_BLK_SIZE ; i++){
	B[i+j*L2_BLK_SIZE] = 
	  A[(submatrix_row*L2_BLK_SIZE) + i + (submatrix_col*L2_BLK_SIZE+j)*h];
      }
    }
  } else if((submatrix_row+1)*L2_BLK_SIZE > h && 
	    (submatrix_col+1)*L2_BLK_SIZE <= w){
    if(verbose) printf("\n=copyL2Block= case 2");
    int diff = (submatrix_row+1)*L2_BLK_SIZE - h;
    int new_rowguard = (L2_BLK_SIZE - diff);
    for(int j = 0; j < L2_BLK_SIZE ; j++){
      for(int i = 0; i < new_rowguard ; i++){
	B[i+j*L2_BLK_SIZE] = 
	  A[(submatrix_row*L2_BLK_SIZE) + i + (submatrix_col*L2_BLK_SIZE+j)*h];
      }
    }
    //fill with zeros for rest
    for(int j = 0; j < L2_BLK_SIZE; j++){
      for(int i = new_rowguard; i < L2_BLK_SIZE; i++){
	B[i+j*L2_BLK_SIZE] = 0.0;
      }
    } 
  } else if((submatrix_row+1)*L2_BLK_SIZE <= h 
	    && (submatrix_col+1)*L2_BLK_SIZE > w){
    if(verbose) printf("\n=copyL2Block= case 3");
    int diff = (submatrix_col+1)*L2_BLK_SIZE - w;
    int new_colguard = (L2_BLK_SIZE - diff);
    for(int j = 0; j < new_colguard ; j++){
      for(int i = 0; i < L2_BLK_SIZE; i++){
       	B[i+j*L2_BLK_SIZE] = 
	  A[(submatrix_row*L2_BLK_SIZE) + i + (submatrix_col*L2_BLK_SIZE+j)*h];
      }
    } 
    //fill with zeros for rest
    for(int j = new_colguard; j < L2_BLK_SIZE; j++){
      for(int i = 0; i < L2_BLK_SIZE; i++){
	B[i+j*L2_BLK_SIZE] = 0.0;
      }
    }
  } else if((submatrix_row+1)*L2_BLK_SIZE > h 
	    && (submatrix_col+1)*L2_BLK_SIZE > w){
    if(verbose) printf("\n=copyL2Block= case 4");
    int diff1 = (submatrix_row+1)*L2_BLK_SIZE - h;
    int new_rowguard = (L2_BLK_SIZE - diff1);
    int diff2 = (submatrix_col+1)*L2_BLK_SIZE - w;
    int new_colguard = (L2_BLK_SIZE - diff2);
    for(int j = 0; j < new_colguard; j++){
      for(int i = 0; i < new_rowguard; i++){
	B[i+j*L2_BLK_SIZE] = 
	  A[(submatrix_row*L2_BLK_SIZE) + i + (submatrix_col*L2_BLK_SIZE+j)*h];
      }
    } 
    //fill with zeros for rest
    for(int j = new_colguard; j < L2_BLK_SIZE; j++){
      for(int i = new_rowguard; i < L2_BLK_SIZE; i++){
	B[i+j*L2_BLK_SIZE] = 0.0;
      }
    }
  }
}
//copy a NOT-YET-transposed L1 submatrix block of L2-sized A to L1-block-sized B
void copyL1Block(const double * const A, const int i, const int j, double * B){
  int verbose = 0;
  if(verbose) printf("=CopyL1Block=");
  for(int j2 = 0; j2 < L1_BLK_SIZE ; j2++){
    for(int i2 = 0; i2 < L1_BLK_SIZE ; i2++){
      B[i2 + j2 * L1_BLK_SIZE] 
	= A[i*L1_BLK_SIZE + (j*L1_BLK_SIZE)*L2_BLK_SIZE + i2 + j2*L2_BLK_SIZE];
    }
  }
}
//given an L2 sized submatrix A, transpose it, one L1 submatrix at a time, into B
void transposeL2(const double * const A, double * B){
  static __attribute__ ((aligned(16))) double b_blockL1[L1_BLK_SIZE*L1_BLK_SIZE];
  static __attribute__ ((aligned(16))) double c_blockL1[L1_BLK_SIZE*L1_BLK_SIZE]; 
  int verbose = 0;
  if(verbose) printf("\n=transposeL2=");
  for(int j = 0; j < MULTIPLE ; j++){
    for(int i = 0; i < MULTIPLE ; i++){
      copyL1Block(A, i ,j , b_blockL1);
      if(verbose){
	printf("\nb_blockL1:\n");
	prettyPrint(b_blockL1, L1_BLK_SIZE, L1_BLK_SIZE);
      }
      transposeL1Size(b_blockL1, c_blockL1);
      if(verbose){
	printf("\nc_blockL1:\n");
	prettyPrint(c_blockL1, L1_BLK_SIZE, L1_BLK_SIZE);
      }
      copyTransposedL1Block(c_blockL1, i, j, B);
      //      printf("\nUpdated L2:\n");
      //prettyPrint(B, L1_BLK_SIZE, L1_BLK_SIZE);
    }
  }
}
//copy a tranposed L1-block sized block B into A which is L2 sized, while 
//transposing submatrix itself amongst the other submatrices
void copyTransposedL1Block(const double * const B, const int submatrix_row, 
			   const int submatrix_col, double * A){
  for(int j = 0; j < L1_BLK_SIZE; j++){
    for(int i = 0; i < L1_BLK_SIZE; i++){
      A[submatrix_col*L1_BLK_SIZE + submatrix_row*L1_BLK_SIZE*L2_BLK_SIZE 
	+ i + j*L2_BLK_SIZE] = 
	B[i+j*L1_BLK_SIZE];
    }
  }
}
//fill another supplied matrix with the tranpose of this one; for fixed size
//only so no need to worry about corner cases
void transposeL1Size(double * const A, double * B){
  for(int j = 0; j < L1_BLK_SIZE; j++){
    for(int i = 0; i < L1_BLK_SIZE; i++){
      B[i+j*L1_BLK_SIZE] = A[j+i*L1_BLK_SIZE];
    }
  }
}
