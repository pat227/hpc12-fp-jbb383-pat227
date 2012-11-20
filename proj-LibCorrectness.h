#ifndef MATRIXCHECKS
#define MATRIXCHECKS
//use these as boolean functions
/* Check for constraint: Q * Qtranspose = I 
   Apply to Q of the A=QR decomposition. Note that Q must be square. I've 
   seen some decompositions called "QR" where Q is not square...these are
   not supported. This function does not store the product matrix, but 
   only computes individual elements at a time and compares against 
   identity matrix.
   m -> # of rows and cols of Q
*/
int IsQbyQtransposeIdentity(const double * const Q, const uint32_t m);

/* Check for constraint: A = QR given QR; need additional args for m&n of Q and 
   R each 
   Function fails fast if any error occurs and doesn't waste time computing rest 
   of the matrix. Function does not store the result of the multiplication, but 
   only compares individually computed elements to elements in A. To support 
   non-square R require more args to specify the dimensions of Q and R each.
   Perhaps structs with bounds would be better and they would know how to multiply
   against each other... 
   Q -> pointer to array of doubles in column-major order that represents Q
   R -> ditto that represents R
   A -> ditto that represents A
   Qm, Qn -> the m x n size of Q
   Rm, Rn -> the m x n size of R
*/
int IsQRequalToA(const double * const Q, const double * const R, 
		 const double * const A, const uint32_t Qm, const uint32_t Qn, 
		 const uint32_t Rm, const uint32_t Rn);

/* Check for constraint: R is upper triangular */
int isUpperTriangular(const double * const M, const uint32_t m, 
		      const uint32_t n);
void print_matrix(const double matrix[], const int m, const int n);
/*not needed here except for testing...might have to use some conditional 
  compilation*/
double identity(const uint32_t m, const uint32_t n);
#endif

