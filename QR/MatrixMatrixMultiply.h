/* Header for Matrix Matrix Multiply */

/* Matrix Matrix Multiply */
int MatrixMatrixMultiply( double *, int, int, double *, int, int, double *); 

/* Code to load and unload blocks of matrix */
void BlockMatrix( const double *, double *, int, int, int, int, int);
void UnBlockMatrix( double *, const double *, int, int, int, int, int);
void CleanMatrix( double *, int, int);


