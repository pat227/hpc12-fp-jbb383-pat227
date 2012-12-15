/* HPC 2012 Final Project: Jacqueline Bush, Paul Torres */

/* Header File: Code in this file is used to check the more complicated code used later. */

/* Code to test Q for orthogonality */
void testOrthogonal( double *, double *, int );

/* Code to test if R is upper triangular */
void testUpperTriangular( double  *, int,int);

/* Code to test Blocked Matrix Matrix Multiply */
void testMatrixMultiply( const double *, const double *, double *, int, int, int, int);

/* Code to test Blocked Matrix Transpose */
void testMatrixTranspose(const double *, int , int , const double *);

//Given Q, R & A, test QR=A
int IsQRequalToA(const double * const Q, const double * const R, 
		 const double * const A, const int Qm, const int Qn, 
		 const int Rm, const int Rn);
