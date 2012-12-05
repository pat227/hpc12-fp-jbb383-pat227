/* HPC 2012 Final Project: Jacqueline Bush, Paul Torres */
/* Header File */

int MatrixMatrixMultiply( double *, int, int, double *, int, int, double *); 
void MatrixTranspose(const double * const A, const int w, const int h, double * B);
int WY(double *, int, int,  double *);
//I NEED THIS in 2 files, yet we have no implementation file of the same name (QR.c)...this needs reorganizing!
/*
void prettyPrint(const double * const A, const int w, const int h);
  printf("\n");
  for(int j = 0; j < w ; j++){
    for(int i = 0; i < h; i++){
      printf("%5.5f ", A[i+j*h]);
    }
    printf("\n");
  }
  printf("\n");
}
*/
