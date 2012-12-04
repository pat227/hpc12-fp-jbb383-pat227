/* HPC 2012 Project  : Jacqueline Bush, Paul Torres */

/* WY Code: 
Given a Matrix A this code outputs the transpose of the orthogonal matrix Q.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "QR.h"
#include <math.h>

/* Code to perform Householder */
void CalculateV( double *, int, int, int, double *);


/*=========================================================================*/

int WY( double *A, int w, int h, double *Qtranspose){

/* Calculate v1 */
double *v = malloc(h);
CalculateV( A, w, h, 0, v); 

/* Allocate W and Y */ 
double *W = malloc (h*h);
double *Y = malloc (h*h);

/* Allocate extra arrays needed. */
double *Ytranspose = malloc(h*h);
double *temp = malloc(h*h);
double *z = malloc( h);
double *Q = malloc(h *h );

/* Fill in First column of W and Y, everything becomes 0. */
for(int i=0; i<h; i++){
	Y[i] = v[i]; 
	W[i] = -2*v[i];
	Ytranspose[i] = 0;
	temp[i] = 0;
	Q[i] = 0;
for(int j=1; j<h; j++){
	Y[i + j*h] =0;
	W[i + j*h] = 0;
	Ytranspose[i + j*h] = 0;
	temp[i + j*h] = 0;
	Q[i + j*h] = 0;
} }


/* Enter Loop */
for( int k = 1; k<h; k++){
	/* Calculate new v */
	CalculateV(A, w, h, k, v);
	
	/* Matrix Transpose Y */
	MatrixTranspose(Y, w, h,Ytranspose);

	/* Matrix Multiply W Y^T */
	MatrixMatrixMultiply(W, h,h, Ytranspose, h, h, temp); 
	
	/* Add I to temp */
	for(int i =0; i<h; i++){
		temp[i + h*i] += 1 ;
	}
	
	/* Scalar Multiply temp by -2 */ 
	for(int i=0; i < h; i++){
	    for(int j=0; j < h; j++){
	      temp[i + j*h] *= -2;
	    }
	  }


	/* Matrix Matrix Multiply temp and v */
	MatrixMatrixMultiply( temp, h, h, v, 1, h, z);


	/* Insert new columns into Y and W */
	for(int i=0; i<h; i++){
	Y[i+ k*w] = v[i]; 
	W[i+ k*w] = z[i];
	}
}

/* Perform Matrix Matrix Multiplication on W Y^T */
	MatrixMatrixMultiply(W, h,h, Ytranspose, h, h, Q); 
	
	/* Add I to get real Q */
	for(int i =0; i<h; i++){
		Q[i + h*i] += 1 ;
	}

/* Perform matrix transpose on Q */
  MatrixTranspose(Q, w, h, Qtranspose);


return 0;
}


void CalculateV( double *A, int w, int h, int coli, double *v){

/* Extract v vector */
for(int i =0; i< coli ; i++){
	v[i] = 0;
	}

for(int i=coli; i<h; i++){
	v[i] = A[(coli + i) + coli*w]; 
}

/*Calculate sign of x*/
int sign =1 ; 

if ( v[1] <0)
	sign = -1;

/*Calculate norm */
double normx=0;
double temp =0;

for( int i =0; i<(h-coli) ; i++){
   temp = v[i];
   temp *= temp;
   normx += temp;
}
   normx = sqrt(normx);

/*Update v */
v[0] += sign*normx;  

/*Calculate new norm */
double normv=0;
temp = 0;
for( int i =0; i<(h-coli) ; i++){
   temp = v[i];
   temp *= temp;
   normv += temp;
}
   normv = sqrt(normv);


/*Update v */
for(int i = 0; i< h-coli; i++)
	v[i] /= normv;

}




