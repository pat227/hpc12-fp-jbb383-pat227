/* HPC 2012 Final Project: Jacqueline Bush, Paul Torres */

/* Header File: Code in this file is used to check the more complicated code used later. */

/* Code to test Q for orthogonality */
void testOrthogonal( double *, double *, int );

/* Code to test if R is upper triangular */
void testUpperTriangular( double  *, int,int);

/* Code to test Blocked Matrix Matrix Multiply */
//this is now being used--not for testing--by QR code, should be moved to utilities
//void dgemm_simple( const double *, const int, const int, const double *, const int, const int, double *);

void testMatrixMultiply( const double *, const double *, double *, int, int, int, int);

/* Code to test Blocked Matrix Transpose */
//ditto...move to utilities
//void simple_transpose(const double *, int , int , double *);
void testMatrixTranspose(const double *, int , int , const double *);
//Given Q, R & A, test QR=A
int IsQRequalToA(const double * const Q, const double * const R, 
		 const double * const A, const int Qm, const int Qn, 
		 const int Rm, const int Rn);
