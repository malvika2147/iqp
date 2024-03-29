/* Library for basic matrix operations over F2 
*  Main operation is to sample a random vector from the null space of M
*  Requires user defined "coin" which can generate random bits.
*/

#include <vector>
#include <functional>

typedef int bit;
typedef std::vector<std::vector<bit>> Matrix;
typedef std::vector<bit> Vector;

void cnot(const Vector &a, Vector &b);
Matrix rref(const Matrix &A);
int rnk(const Matrix &A);

// Return random vector from null space
Vector sample_from_null(const Matrix &A, std::function<bit()> &coin);
// Return a random solution to the system of linear equations Ax=y
Vector sample_from_solutions(const Matrix &A, const Vector &y, std::function<bit()> &coin);

// Get all vectors from null space
Matrix get_entire_null(const Matrix &A);

// Compute Mv 
Vector mul(const Matrix &M, const Vector &v);
// Matrix multiplication
Matrix mul(const Matrix &A, const Matrix &B);

// Compute dot product
bit dot(const Vector &a, const Vector &b); 
// Check if 0 vector
bool check_zero(const Vector &a);

Matrix transpose(const Matrix &A);

// Returns random mxn matrix
Matrix gen_random_matrix(int m, int n, std::function<bit()> &coin);
Vector gen_random_vector(int n, std::function<bit()> &coin);


void print(const Matrix &M);
// Default orientation is row vector
void print(const Vector &v, bool orient=0);
