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


/*=========================================================================*/

int WY( double *A, int h, int w, double *Q){

/* Current Householder vector */
double *v = malloc( h * sizeof(double));
double *z = malloc( h * sizeof(double));

/* Temporary Arrays */
double *temp = malloc(h * h* sizeof(double));
double *temp1 = malloc(h * w * sizeof(double));
double *temp2 = malloc(h * h * sizeof(double));
double *temp3 = malloc(h* h * sizeof(double));
double *temp4 = malloc( h * h * sizeof(double));
double *temp5 = malloc( h* h * sizeof(double));
double *temp6 = malloc(h * h * sizeof(double));
double *temp7 = malloc(h * h * sizeof(double));

/* W and Yt Arrays. Note: Yt = Y^t */
double *W = malloc( h* w* sizeof(double));
double *Yt = malloc( w*h * sizeof(double)); 



/*============================  First iteration ==========================================*/


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




MatrixMatrixMultiply( v, h, 1, v, 1, h , temp2);
/* Scalar multiplication by -2 */
	for(int i=0; i< h; i++){
	for( int j=0; j<h; j++){
		temp2[i + j*h] *=-2 ;
	} }
for(int i=0; i<h; i++){
	temp2[ i + i*h] += 1;
}


	CalculateQ( W, Yt, h, w, temp);

	MatrixMatrixMultiply( temp, h,h,A, h, w, temp1);
	
	printf(" R k %d : \n", 0); 
	prettyPrint(temp1, h, w);


/*============================= Second Iteration ==========================================*/

	CalculateV( temp1, h, w, 1, v);
	prettyPrint(v, h, 1);

MatrixMatrixMultiply( v, h, 1, v, 1, h, temp3);
/* Scalar multiplication by -2 */
	for(int i=0; i< h; i++){
	for( int j=0; j<h; j++){
		temp3[i + j*h] *=-2 ;
	} }

for(int i=0; i<h; i++){
	temp3[ i + i*h] += 1;
}



MatrixMatrixMultiply(temp2, h,h, temp3, h,h, temp4);

	/* Scalar multiplication by -2 */
	for(int i=0; i< h; i++){
	for( int j=0; j<h; j++){
		temp[i + j*h] *=-2 ;
	} }

	MatrixMatrixMultiply(temp, h, h, v, h ,1, z);
	prettyPrint(z, h, 1); 

	/* Fill in k row or column of W and Yt */
	for(int i=0; i< h; i++){
	W[i+ h] = z[i];
	Yt[1+ i*w] = v[i];
	}

	CleanMatrix(temp, h, h);

	CalculateQ( W, Yt, h, w, temp);
	prettyPrint(temp, h, h);
	prettyPrint(temp4, h, h);

  	CleanMatrix(temp1, h,w );      

	MatrixTranspose( temp, h, h, temp5);
	
	MatrixMatrixMultiply( temp5, h,h,A, h, w, temp1);
	
	printf(" R k %d : \n", 1); 
	prettyPrint(temp1, h, w);

/*=============================  Thrid iteration ======================================== */

	CalculateV( temp1, h, w, 2, v);
	printf("v = \n");
	prettyPrint(v, h, 1);

MatrixMatrixMultiply( v, h, 1, v, 1, h, temp6);
/* Scalar multiplication by -2 */
	for(int i=0; i< h; i++){
	for( int j=0; j<h; j++){
		temp6[i + j*h] *=-2 ;
	} }

for(int i=0; i<h; i++){
	temp6[ i + i*h] += 1;
}



      MatrixMatrixMultiply(temp4, h,h, temp6, h,h, temp7);

	/* Scalar multiplication by -2 */
	for(int i=0; i< h; i++){
	for( int j=0; j<h; j++){
		temp[i + j*h] *=-2 ;
	} }

	MatrixMatrixMultiply(temp, h, h, v, h ,1, z);
	prettyPrint(z, h, 1); 

	/* Fill in k row or column of W and Yt */
	for(int i=0; i< h; i++){
	W[i+ h] = z[i];
	Yt[1+ i*w] = v[i];
	}

	CleanMatrix(temp, h, h);

	CalculateQ( W, Yt, h, w, temp);
	
        prettyPrint(temp, h, h);
	prettyPrint(temp7, h, h);

  	CleanMatrix(temp1, h,w );      
	


	MatrixTranspose( temp7, h, h, temp5);
	
	MatrixMatrixMultiply( temp5, h,h,A, h, w, temp1);
	
	printf(" R k %d : \n", 1); 
	prettyPrint(temp1, h, w);






return 0;
}



void CalculateQ( double *W, double *Yt, int h, int w, double *temp){

/* Calculates temp =  I + W Y^T */

MatrixMatrixMultiply( W, h, w, Yt, w, h, temp);

for(int i=0; i<h; i++){
	temp[ i + i*h] += 1;
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






