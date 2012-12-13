#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "Utilities.h"
#include "test.h"

int MatrixVectorMultiply(const double * const A, const int hA, const int wA, 
			 const double * const B, double *C);
