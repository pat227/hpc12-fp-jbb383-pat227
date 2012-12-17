/* HPC 2012 Project  : Jacqueline Bush, Paul Torres */

/* WY Code: 
Given a Matrix A this code outputs the transpose of the orthogonal matrix Q.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
//#include "MatrixMatrixMultiply.h"
//#include "MatrixTranspose.h"
#include "WY.h"
#include "test.h"
#include "Utilities.h"
#include <math.h>


/*=========================================================================*/

int WY( double *A, int h, int w, double *Q, double *Qt, double *R){

/* Initialize Matrices and vectors */
  double *v = malloc( h * sizeof(double));
  double *z = malloc( h * sizeof(double));
	double *a1 = malloc( h * sizeof(double));
	double *a2 = malloc( h * sizeof(double)); 
	double *Yt = malloc( w * h * sizeof(double));
	double *W = malloc( h * w * sizeof(double));


/* Set Yt and W to zero matrix */
	CleanMatrix(Yt, w, h);
	CleanMatrix(W, h, w);

/*===================== Initial steps =========================== */

/* Copy a_0 into a1 */
	for(int i=0; i<h; i++){
		a1[i] = A[i]; 
	}

/* Calculate initial v */
	CalculateV( a1, h, 0, v);	

/* Fill in first row or column */	
	for(int i = 0; i < h; i++){
		W[i] = -2*v[i];
		Yt[i*w] = v[i];
	}

/* Calculate initial Q = I + W Yt */
	CalculateQ(W, Yt, h, w, Q);
		
/* Calculate Q^T */
	simple_transpose(Q, h, h, Qt);
	//MatrixTranspose(Q, h, h, Qt);

/* Set n equal to min(w, h). */
	int n = w;
	if( w > h)
		n= h;
	
/* ========================== Enter Loop ========================= */
for(int k=1; k<n; k++){


	/* Copy a_k relavant entries into a1 and make such it is not a zero vector*/
	double check= 0;

	for(int i=0; i<h; i++){
		a1[i] = A[i+k*h];
		if( i > k) 
		check += A[i +k*h];
	}

	/* Don't continue in loop if we have reached a zero column */
	if ( check == 0 ){
		break;
	}

	/* Update a_k */
	//MatrixMatrixMultiply(Qt, h, h, a1, h, 1, a2);
	dgemm_simple(Qt, h, h, a1, h, 1, a2);

	/* Calculate kth v */
	CalculateV( a2, h, k, v);	

	/* Calculate kth z */
	//MatrixMatrixMultiply(Q,h,h,v, h,1, z);
	dgemm_simple(Q, h, h, v, h, 1, z);	

	/* Fill in the kth row or column */	
	for(int i = 0; i < h; i++){
		W[i+ k*h] = -2*z[i];
		Yt[k +i*w] = v[i];
	}

	/* Calculate kth Q = I + W Yt */
	CalculateQ(W, Yt, h, w, Q);
	
	/* Calculate Q^T */
	simple_transpose(Q, h, h, Qt);
	//MatrixTranspose(Q, h, h, Qt);

}	

/*========================= Calculate R ==============================*/
//MatrixMatrixMultiply(Qt, h,h,A , h ,w , R);
dgemm_simple(Qt, h, h, A, h, w, R);

/*====================== Test Code =================================*/

/* Uncomment if you want to test this code */
//testOrthogonal(Q, Qt, h);
//testUpperTriangular(R, h, w);



/*======================== Clean Up  ====================================*/
//double free error if we try next 2 lines - BUG someplace
//free(v); 
//free(z); 
free(a1); 
free(a2); 
free(W); 
free(Yt);

return 0;
}



void CalculateQ( double *W, double *Yt, int h, int w, double *Q){

/* Calculates temp =  I + W Y^T */

dgemm_simple(W, h, w, Yt, w, h, Q);
//MatrixMatrixMultiply( W, h, w, Yt, w, h, Q);

for(int i=0; i<h; i++){
	Q[ i + i*h] += 1;
}

}

void CalculateV( double *A, int h, int coli, double *v){
/* Calculates desired housholder vector */

/* Extract v vector */
	for(int i =0; i< coli ; i++){
		v[i] = 0;
	}

	for(int i=coli; i<h; i++){
		v[i] = A[i]; 
	}

/*Calculate sign of x*/
	int sign =1 ; 

	if ( v[coli] <0)
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
