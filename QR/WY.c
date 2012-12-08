/* HPC 2012 Project  : Jacqueline Bush, Paul Torres */

/* WY Code: 
Given a Matrix A this code outputs the transpose of the orthogonal matrix Q.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "MatrixMatrixMultiply.h"
#include "MatrixTranspose.h"
#include "WY.h"
#include "test.h"
#include "Utilities.h"
#include <math.h>

/* Currently works for k = 2, We need m greater than or equal to n */


/*=========================================================================*/

int WY( double *A, int h, int w, double *Q){

/* Current Householder vector */
double *v = malloc( h * sizeof(double));

/* Current z vector */
double *z = malloc( h * sizeof(double));

/* Qt Array */
double *Qt = malloc(h * h * sizeof(double));


/* R array */
double *R = malloc( h* w * sizeof(double));

/* W and Yt Arrays. Note: Yt = Y^t */
double *W = malloc( h* w* sizeof(double));
double *Yt = malloc( w*h * sizeof(double)); 

/* Initial Householder vector*/
CalculateV( A, h, w, 0, v);

/* Clean W and Yt */
CleanMatrix(W, h, w);
CleanMatrix(Yt, w, h);

/* Fill in first row or column of W and Yt */
for(int i=0; i< h; i++){
	W[i] = -2 * v[i];
	Yt[i*w] = v[i];
}

/* Calculate Temp */	
	CleanMatrix(Q, h, h);
	CalculateQ( W, Yt, h, w, Q);

		
/* Calculate new R */
	MatrixTranspose( Q, h,h, Qt);	
	MatrixMatrixMultiply( Qt, h,h,A, h, w, R);
	
	printf(" R k=%d : \n", 0); 
	prettyPrint(R, h, w);


/*********************** Start Loop *********************************/

for( int k= 1; k< w; k++){

	CalculateV( R, h, w, k, v);

/* Calculate new z */

	/* Scalar multiplication by -2 */
	for(int i=0; i< h; i++){
	for( int j=0; j<h; j++){
		Q[i + j*h] *=-2 ;
	} }

	MatrixMatrixMultiply(Q, h, h, v, h ,1, z); 

/* Fill in k row or column of W and Yt */
	for(int i=0; i< h; i++){
	W[i+ k*h] = z[i];
	Yt[k+ i*w] = v[i];
	}	

/* Calculate Temp */	
	CleanMatrix(Q, h, h);
	CalculateQ( W, Yt, h, w, Q);

		
/* Calculate new R */
	CleanMatrix(R, h, w);
	MatrixTranspose( Q, h,h, Qt);	
	MatrixMatrixMultiply( Qt, h,h,A, h, w, R);
	
	printf(" R=k %d : \n", k); 
	prettyPrint(R, h, w);


}

}



void CalculateQ( double *W, double *Yt, int h, int w, double *Q){

/* Calculates Q =  I + W Y^T */

MatrixMatrixMultiply( W, h, w, Yt, w, h, Q);

for(int i=0; i<h; i++){
	Q[ i + i*h] += 1;
}

}

void CalculateV( double *A, int h, int w, int coli, double *v){
/* Calculates desired housholder vector */

/* Extract v vector */
	for(int i =0; i< coli ; i++){
		v[i] = 0;
	}

	for(int i=coli; i<h; i++){
		v[i] = A[(i) + coli*h]; 
	}

/*Calculate sign of x*/
	int sign =1 ; 

	if ( v[0] <0)
		sign = -1;

/*Calculate norm */
	double normx=0;
	double temp =0;

	for( int i =0; i< h ; i++){
   		temp = v[i];
   		temp *= temp;
   		normx += temp;
	}
   	normx = sqrt(normx);

/*Update v */
	v[coli] += sign*normx;  

/*Calculate new norm */
	double normv=0;
	temp = 0;
	for( int i =0; i< h ; i++){
   		temp = v[i];
   		temp *= temp;
   		normv += temp;
	}
   	normv = sqrt(normv);

/*Update v */
	for(int i = 0; i< h; i++)
		v[i] /= normv;


}




