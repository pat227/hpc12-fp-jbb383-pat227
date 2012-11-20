#ifndef MATRIXCHECKS
#define MATRIXCHECKS
//use these as boolean functions
/* Check for constraint: A * A transpose = I */
int IsAbyAtransposeIdentity(const double * const A, const uint32_t m, 
			    const uint32_t n);

/* Check for constraint: A = QR given QR; need additional args for m&n of Q and 
   R each */
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

