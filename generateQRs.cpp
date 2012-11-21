#include<iostream>
#include<armadillo> //see http://arma.sourceforge.net/, better than lapack alone
using namespace std;
using namespace arma;

void print_matrixWithMorePrecision(const mat &A);

int main(int argc, char** argv){
  /*
  mat A = randu<mat>(4,4);
  mat B = randu<mat>(4,4);
  cout << "A is: " << endl << A;
  cout << "B is: " << endl << B << endl;
  cout << "A*B is: " << endl;
  cout << A*B << endl;
  cout << "B transposed is:" << endl;
  cout << B.t() << endl;
  cout << "A * Btransposed is:" << endl;
  cout << A*B.t() << endl;
    */

    /*
      4 5 4
      4 8 8
      2 6 6
      1 7 4
      6 1 7
     */
  mat A = mat(6,3);
  A(0,0) = 4;
  A(0,1) = 5;
  A(0,2) = 4;
  A(1,0) = 4;
  A(1,1) = 8;
  A(1,2) = 8;
  A(2,0) = 2;
  A(2,1) = 6;
  A(2,2) = 6;
  A(3,0) = 1;
  A(3,1) = 7;
  A(3,2) = 4;
  A(4,0) = 6;
  A(4,1) = 1;
  A(4,2) = 7;
  A(5,0) = 9;
  A(5,1) = 0;
  A(5,2) = 2;
  cout << "Matrix A:" << endl << A << endl;
  cout << "QR factorization according to Armadillo++ lib (using lapack under the hood):" << endl;
  mat Q, R;
  qr(Q,R,A);
  cout << "Q:" << endl;
  //does not have desired effect
  //cout.precision(9);
  //cout << Q;  
  //cout << "R:" << endl << R;
  //mat * const pMatQ = &Q;
  //const mat * const pMatR = &R;
  print_matrixWithMorePrecision(Q);
  cout << "R:" << endl;    
  print_matrixWithMorePrecision(R);
  return 0;
}


void print_matrixWithMorePrecision(const mat &A){
  cout.precision(9);
  cout.width(9);
  double temp = 0.0;
  for(int i = 0; i < A.n_rows; i++){
    for(int j = 0; j < A.n_cols; j++){
      temp = A(i,j);
      cout<<(temp);
      cout<<("\t");
    }
    cout<<("\n");
  }
}
