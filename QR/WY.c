/* HPC 2012 Project  : Jacqueline Bush, Paul Torres */

/* WY Code: 
Given a Matrix A this code outputs the transpose of the orthogonal matrix Q.
*/


/* Note: Not quite working yet.... */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "MatrixMatrixMultiply.h"
#include "MatrixTranspose.h"
#include "WY.h"
#include "test.h"
#include "Utilities.h"
#include <math.h>


/*=========================================================================*/

int WY( double *A, int h, int w, double *Q){

/* Allocate vectors v and z */
double *v = malloc(h*sizeof(double));
double *z = malloc(h*sizeof(double));

/* Allocate matrices Y and W, We will store Y in transposed form to save on computational costs*/
double *Y = malloc(h*h*sizeof(double));
double *W = malloc(h*h*sizeof(double));
double *Qt = malloc(h*h*sizeof(double));


/* Allocate temp matrices */
double *temp = malloc(h*w*sizeof(double));

/* Initilize Arrays to zero */
CleanMatrix(v, h, 1);
CleanMatrix(z,h, 1);
CleanMatrix(W, h,h);
CleanMatrix(Y, h,h);
CleanMatrix(temp,h,w);
CleanMatrix(Qt, h,h);

prettyPrint(A,h,w);


/* Initial v and z */
CalculateV( A, h, w, 0, v);

/* Fill in 1st column of W and 1st row of Y */
	for(int i = 0; i<h; i++){
		W[i] = -2*v[i];
		Y[i*h] = v[i];
	}

    printf(" W = \n");
	prettyPrint( W, h,h);
	printf("Y = \n");
	prettyPrint(Y, h,h);

/* Clean Vectors v and z */
	CleanMatrix(v, h, 1);
	CleanMatrix(z, h, 1);

for(int k =1; k<1; k++){
	/* Calculate v */
	CalculateV( A, h, w, k, v); 

	/* Calculate z */
	CalculateZ(W, Y, v, h, z );	
	
	/* Fill in kth column of W and kth row of Y */
	for(int i = 0; i<h; i++){
		W[i+ k*h] = z[i];
		Y[k + i*h] = v[i];
	}
	
   printf(" W = \n");
	prettyPrint( W, h,h);
	printf("Y = \n");
	prettyPrint(Y, h,h);


	/* Clean Vectors v and z */
	CleanMatrix(v, h, 1);
	CleanMatrix(z, h, 1);
}

MatrixMatrixMultiply(W, h,h, Y, h,h, Q);


prettyPrint(Q, h,h);

for(int i =0; i< h; i++){
		Q[i + i*h] += 1;
	}

printf("Q = \n");
prettyPrint(Q, h,h);
simple_transpose( Q, h, h, Qt );
printf("Qt = \n");
prettyPrint(Qt, h,h);
 
 MatrixMatrixMultiply(Qt, h,h, A, h, w, temp);

printf("R = \n");
prettyPrint(temp, h, w);

return 0;
}


void CalculateZ( double *W, double *Y, double *v, int h, double *z){
/* Calculates z = -2 ( I + WY) v */

/* Allocate and clean temporay array */
	double *temp = malloc(h*h*sizeof(double));
	CleanMatrix(temp, h,h);

/*  W*Y (Note Y = Y^T theoretically)*/

   MatrixMatrixMultiply( W, h,h, Y, h,h, temp);

printf(" WY = \n");
prettyPrint(temp, h,h);

/* I + W*Y */
	for(int i =0; i< h; i++){
		temp[i + i*h] += 1;
	}
printf(" temp = \n");
prettyPrint(temp, h, h);

/* (I + W*Y)*v */ 
	MatrixMatrixMultiply( temp, h, h, v, h, 1 , z);

/* -2 ( I + W*Y)*v */
	for(int i =0; i< h; i++){
		z[i] *= -2;
	}	

}



void CalculateV( double *A, int h, int w, int coli, double *v){

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




