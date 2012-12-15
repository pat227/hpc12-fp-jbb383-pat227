/* HPC 2012 Final Project: Jacqueline Bush, Paul Torres */
/* Header file for random functions */

void prettyPrint(const double * const A, const int m, const int n);

/* Code to load and unload blocks of matrix */
void BlockMatrix( const double *, double *, int, int, int, int, int);
void UnBlockMatrix( double *, const double *, int, int, int, int, int);
void CleanMatrix( double *, int, int);
void BlockVectorForMatrixVector(const double *inA, double *outA, const int hA, 
				const int b, const int i_bloc);
void UnBlockVectorForMatrixVector(double *outA, const double *inA, const int hA,
				  const int b, const int i_bloc);
void writetofile(const char * const fname, int m, int n, int iterations, double dependent);
void writetofile2(const char * const fname, int m, int n, double dependent);
//imports from test.h b/c no longer being used merely for testing but for code
void dgemm_simple( const double *, const int, const int, const double *, const int, const int, double *);
void simple_transpose(const double *, int , int , double *);
