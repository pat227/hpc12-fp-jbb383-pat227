/* HPC 2012 Project: Jacqueline Bush, Paul Torres */

/* Header file for Matrix Transpose */


void copyTransposedL2Block(const double * const B, const int submatrix_row , 
			   const int submatrix_col, double * A, const int w, const int h);
void copyL2Block(const double * const A, const int submatrix_row , 
		 const int submatrix_col, double * B, const int w, const int h);
void copyL1Block(const double * const A, const int i, const int j, double * B);
void transposeL2(const double * const A, double * B);
void copyTransposedL1Block(const double * const B, const int submatrix_row, 
			   const int submatrix_col, double * A);
void transposeL1Size(double * const A, double * B);
void MatrixTranspose(const double * const A, const int h, const int w, double * B);
