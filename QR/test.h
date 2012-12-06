/* HPC 2012 Final Project: Jacqueline Bush, Paul Torres */

/* Header File: Code in this file is used to check the more complicated code used later. */


/* Code to test Blocked Matrix Matrix Multiply */
void dgemm_simple( const double *, const int, const int, const double *, const int, const int, double *);
void testMatrixMultiply( const double *, const double *, double *, int, int, int, int);

/* Code to test Blocked Matrix Transpose */
void simple_transpose(const double *, int , int , double *);
void testMatrixTranspose(const double *, int , int , const double *);
