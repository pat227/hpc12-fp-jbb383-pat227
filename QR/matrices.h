
struct matrix{
  int height;
  int width;
  double * elements;
};

//Functions that operate upon matrices
//basic functions
void init(struct matrix * this, int h, int w);
void destroy(struct matrix * this);
void copyMatrix(struct matrix * this, const struct matrix * const other);
double getElement(const struct matrix * const this, const int row, 
			const int col);
void fillWithRandomElements(struct matrix * const m, int upper, int lower);
void prettyPrint(const struct matrix * const this);
int setElement(struct matrix * this, const int row, const int col, 
		     const double value);
void setToIdentity(struct matrix * this);
void swapMatrix(struct matrix * m, struct matrix * other);
void zero(struct matrix * this);
//more advanced functions
void add(const struct matrix * const this, 
	       const struct matrix * const other, struct matrix * answer);
void extractVector(const struct matrix * const this, int row, int col, 
			 struct matrix * vector); 
int matrixMultiply(const struct matrix * const this, 
			const struct matrix * const B, 
			struct matrix * const answer);
void scalarMultiply(struct matrix * this, const double value);
void subtract(const struct matrix * const this, 
		    const struct matrix * const other, struct matrix * answer);
void subtractFromRightBottomMost(struct matrix * const this, 
				const struct matrix * const other);
void transpose(const struct matrix * const this, 
		     struct matrix * transpose);
//compute norms of 1xn or mx1 matrix
double normOfColVector(const struct matrix * const m, const int col, const int rowstart);
double normOfVector(const struct matrix * const m, const int start);
