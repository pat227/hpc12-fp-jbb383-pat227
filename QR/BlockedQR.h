/* HPC 2012 Final Project: Jacqueline Bush and Paul Torres */

/* Header for Blocked QR code */

void BlockedQR( double *A, int h, int w, double *Q, double *R);
void BlockQ(const double *inQ, double *outQ, int hQ, int wQ, int b, int i_bloc1, int j_bloc1, int i_bloc2, int j_bloc2 );
void UnBlockQ( const double *inQ, double  *outQ, int hQ, int wQ, int b, int i_bloc1, int j_bloc1, int i_bloc2, int j_bloc2);
void Block( double *outA, const double *inA, int h, int w, int b, int i_bloc1, int i_bloc2, int j_bloc);
