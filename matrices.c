#include<stdio.h> //fprintf, printf
#include<stdint.h> //exact types
#include<math.h> //sqrt
#include<stdlib.h> //memory, rand
#include<stdarg.h> //var args
#include"matrices.h"
//Computes the sum of 2 matrices into a 3rd provided empty matrix, provided the
//matrices are of equal size
void add(const struct matrix * const m, const struct matrix * const other, 
	 struct matrix * answer){
  if(m->height != other->height || m->width != other->width){
    printf("\nMatrices not of same size for addition; doing nothing...\n");
    return;
  }
  //ensure answer can hold the answer, and if not, fix that
  if(m->height != answer->height || m->width != answer->width){
    destroy(answer);
    init(answer, m->height, m->width);
    answer->width = m->width;
    answer->height = m->height;
  }
  double sum = 0.0;
  for(int i = 0; i < m->height; i++){
    for(int j = 0; j < m->width; j++){
      sum = getElement(m,i,j) + getElement(other,i,j);
      setElement(answer,i,j, sum);
    }
  }
}
//copy element for element
void copyMatrix(struct matrix * m, const struct matrix * const other){
  //ensure they are the same size or else change our size to match
  if(m->height != other->height || m->width != other->width){
    destroy(m);
    init(m, other->height, other->width);
    m->height = other->height;
    m->width = other->width;
  }
  for(int i=0; i < m->height; i++){
    for(int j=0; j < m->width; j++){
      setElement(m,i,j,getElement(other,i,j));
    }
  }
}
//free elements of any matrix
void destroy(struct matrix * m){
  m->width = 0;
  m->height = 0;
  free(m->elements);
}
//Fill another matrix with the col of this matrix, starting with a specified 
//row to the bottom only treat every vector as a matrix of 1 column when we 
//extract it (we never need to extract row vectors)
void extractVector(const struct matrix * const m, const int row, 
		   const int col, struct matrix * vector){
  //ensure we can store all the elements or else resize
  if(vector->width != 1 || vector->height != (m->height-row)){
    destroy(vector);
    init(vector, 1, (m->height-row));
    vector->width = 1;
    vector->height = m->height - row;
  }
  int vecrow = 0;
  for(int i = row; i < m->height; i++){
    setElement(vector, vecrow++, 0, getElement(m,i,col));
  }
}
//fills all elements with random values given an upper and lower bound for those values
void fillWithRandomElements(struct matrix * m, uint32_t upper, uint32_t lower){
  double value = 0.0;
  srand(1);
  for(int j = 0; j < m->width; j++){
    for(int i = 0; i < m->height; i++){
      value = (rand()/(double)RAND_MAX) * (upper - lower) + lower;
      setElement(m, i, j, value);
    }
  }
}

//returns value of an element
double getElement(const struct matrix * const m, const int row, const int col){
  if(row > m->height || col > m->width || row < 0 || col < 0){
    printf("\nOut of bounds while trying to get an element; client wants row %d and col %d but height is %d and width is %d \n", row, col, m->height, m->width);
    return 0;
  }
  return m->elements[row + col*m->height];
}
//initialize a matrix struct for a given set of dimensions, w and h
void init(struct matrix * m, const int w, const int h){
  double * p = malloc(sizeof(double)*w*h);
  if(p==NULL){
    printf("\nCould not allocate memory for a matrix!\n");
    abort();
  }
  m->elements = p;
  m->width = w;
  m->height = h;
}
//given two matrices A & B and a blank answer matrix, fills the 3rd one with the
//product of AB. Returns 0 if not compatible; 1 otherwise.
int matrixMultiply(const struct  matrix * const A, const struct matrix * const B, 
	     struct matrix * const answer){
  //ensure matrices are compatible & so is answer matrix, if latter is not, fix it
  if(A->width != B-> height){
    printf("\nMatrices are not multiplication-compatible; A has width %d, B has height %d.\n", 
	   A->width, B->height);
    return 0;
  }
  if(answer->width != B-> width || answer->height != A->height){
    destroy(answer);
    init(answer,A->height, B->width);
    answer->width = B->width;
    answer->height = A->height;
  }
  double temp = 0.0;
  for(int i=0; i < answer->height; i++){
    for(int j=0; j < answer->width; j++){
      for(int k = 0; k < A->width; k++){
	temp += getElement(A,i,k) * getElement(B,k,j);
      }
      setElement(answer,i,j,temp);
      temp=0.0;
    }
  }
  return 1;
}
//computes norm of next count elements in a column (assuming col-major storage)
//must be given the # of elements in a column
/*
double norm2(const int count, const double * const a, const int width){
  double sum = 0.0;
  for(int i = 0; i < count; i++){
    sum += *(a+width*i);
  }
  return (sum / (double)count);
}
*/
//compute the norm of a col-vector of a submatrix (ie, compute the norm of the 
//elements in col j starting at row i to the bottom of the matrix)
double normOfColVector(const struct matrix * const m, const int col, const int rowstart){
  if(rowstart > m->height || col > m->width || rowstart < 0 || col < 0){
    printf("\nOut of bounds while trying to compute norm...\n");
    return 0;
  } 
  double norm = 0.0;
  double temp = 0.0;
  for(int i = rowstart; i < m->height; i++){
    temp = getElement(m,i,col);
    temp *= temp;
    norm += temp;
  }
  return sqrt(norm);
}
//compute the norm of any single col-sized or single row-sized vector-matrix, 
//starting at any position
double normOfVector(const struct matrix * const m, const int start){
  if(start > m->height || start < 0){
    printf("\nOut of bounds while trying to compute norm...\n");
    return 0;
  } 
  double norm = 0.0;
  double temp = 0.0;
  if(m->height == 1){
    for(int i = start; i < m->width; i++){
      temp = getElement(m,0,i);
      temp *= temp;
      norm += temp;
    }
  } else if(m->width == 1){
    for(int i = start; i < m->height; i++){
      temp = getElement(m,i,0);
      temp *= temp;
      norm += temp;
    }
  } else {
    printf("\nError; norm function invoked upon non-vector; returning 1.0...");
    return 1.0;
  }
  return sqrt(norm);  
}
//print any sized matrix
void prettyPrint(const struct matrix * const m){
  printf("\n");
  for(int i = 0; i < m->height ; i++){
    for(int j = 0; j < m->width ; j++){
      printf("%5.5f ", getElement(m,i,j));
    }
    printf("\n");
  }
  printf("\n");
}
//multiply all values by a real value; division is possible with this method, 
//of course.
void scalarMultiply(struct matrix * m, const double value){
  for(int i=0; i < m->height; i++){
    for(int j=0; j < m->width; j++){
      m->elements[i+j*m->height] *= value;
    }
  }
}
//return 1 if done, 0 if out of bounds
int setElement(struct matrix * m, const int row, const int col, 
	       const double value){
  if(row > m->height || col > m->width || row < 0 || col < 0){
    printf("\nOut of bounds while trying to set an element...\n");
    return 0;
  }
  m->elements[row + col*m->height] = value;
  return 1;
}
//set any sized matrix to identity matrix
void setToIdentity(struct matrix * m){
  for(int i=0; i < m->height; i++){
    for(int j=0; j < m->width; j++){
      if(i == j){
	m->elements[i+j*m->height] = 1.0;} 
      else {
	m->elements[i+j*m->height] = 0.0;
      }
    }
  } 
}
//not really needed but this method saves the time of having to scalar multiply
//an entire matrix by -1 to accomplish subtraction with add()
void subtract(const struct matrix * const m, 
	      const struct matrix * const other, 
	      struct matrix * answer){
  if(m->height != other->height || m->width != other->width){
    printf("\nMatrices not of same size for addition; doing nothing...\n");
    return;
  }
  //ensure answer can hold the answer
  if(m->height != answer->height || m->width != answer->width){
    destroy(answer);
    init(answer, m->height, m->width);
    answer->width = m->width;
    answer->height = m->height;
  }
  double sum = 0.0;
  for(int i = 0; i < m->height; i++){
    for(int j = 0; j < m->width; j++){
      sum = getElement(m,i,j) - getElement(other,i,j);
      setElement(answer,i,j, sum);
    }
  }
}
/*Computes Householder Reflector by subtracting a given matrix from another 
  matrix that is the identity matrix. If smaller, the the other matrix is 
  subtracted from right-most and bottom-most elements only of the identity matrix;
  this function is only of use for Householder transformation computations
  m     -> identity matrix that becomes the Householder Reflector
  other -> the matrix whose values are subtracted from m  
*/
void subtractFromRightBottomMost(struct matrix * const m, 
			  const struct matrix * const other){
  if(m->height < other->height || m->width < other->width){
    printf("\nMatrices of incompatible size for computing H; doing nothing...\n");
    return;
  }
  double sum = 0.0;
  int startrow = m->height - other->height;
  int startcol = m->width - other->width;
  int i2 = 0;
  int j2 = 0;
  for(int i = startrow; i < m->height; i++, i2++){
    for(int j = startcol; j < m->width; j++, j2++){
      sum = getElement(m,i,j) - getElement(other,i2,j2);
      setElement(m,i,j, sum);
    }
    j2 = 0;
  }
}
//fill another supplied matrix with the tranpose of this one
void transpose(const struct  matrix * const m, struct matrix * transpose){
  if(transpose->width != m->height  || transpose->height != m->width){
    destroy(transpose);
    init(transpose, m->width, m->height);
    transpose->width = m->height;
    transpose->height = m->width;
  }
  for(int j = 0; j < m->width; j++){
    for(int i = 0; i < m->height; i++){
      setElement(transpose, j, i, getElement(m,i,j));
    }
  }
}
//zero out all entries in a matrix of any size
void zero(struct matrix * m){
  for (int i=0; i < m->height; i++){
    for(int j=0; j < m->width; j++){
      setElement(m,i,j, 0.0);
    }
  }
}
