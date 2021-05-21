#include "matrix.h"
#include <algorithm>
#include <iostream>
#include <cassert>

using namespace std;

// b += a
void cnot(const Vector &a, Vector &b) {
    int n = a.size();
    for(int i = 0; i < n; ++i) b[i] ^= a[i];
}

bool check_zero(const Vector &a){
    bit ans=0;
    for(int i = 0; i < a.size(); ++i) ans |= a[i];
    return (!ans);
}

Matrix rref(const Matrix &A){
    int n,m;
    m = A.size();
    assert(m);
    n = A[0].size();

    Matrix M = A;

    for(int i = 0,lead=0; i < m; ++i) {
        if(lead >= n) break;
        int k = i;

        // find row with leading 1
        while(M[k][lead] == 0) {
            ++k;
            if(k == m) {
                k=i;
                if(++lead >= n){
                    k=-1;
                    break;
                }
            }
        }

        if(k < 0) break;

        assert(lead < n && M[k][lead]);
        // swap row with leading 1
        for(int j = 0; j < n; ++j) swap(M[i][j],M[k][j]);

        // cnot the remaining rows
        for(k = 0; k < m; ++k) {
            if(k == i) continue;
            if(M[k][lead]) cnot(M[i],M[k]);
        }
    }

    return M;
}

int rnk(const Matrix &A) {
    Matrix M = rref(A);

    int m = M.size();
    assert(m); 
    int n = M[0].size();

    int rank=0;
    for(int i = 0,j=0; i < m; ++i,++j) {
        while(j < n && !M[i][j]) ++j;
        if(j >= n) break;
        ++rank;
    } 

    return rank;
}


vector<bool> get_pivots(const Matrix &M){
    int m = M.size();
    int n = M[0].size();

    vector<bool> lead(n, false);
    for(int i = 0,j=0; i < m; ++i,++j) {
        while(j < n && !M[i][j]) ++j;
        if(j >= n) break;
        lead[j] = 1;
    }
    
    return lead;
}

// takes matrix in rref form and solves system of eq
Vector get_random_solution(const Matrix &M, const Vector &y, std::function<bit()> &coin ){
    int m = M.size();
    int n = M[0].size();

    vector<bool> lead = get_pivots(M);
    Vector sample(n);
    Vector cur(m,0);

    // fix the free variables
    for(int j = 0; j < n; ++j){
        if(lead[j]) {
            continue;
        }
        sample[j] = coin();
        // add sample[j] * column[j] to the answer
        for(int i = 0; i < m; ++i) cur[i] ^= (sample[j]&M[i][j]);
    }
    // calculate the other variables, we want cur to equal to y1
    for(int j =0,k=0 ;j < n; ++j){
        if(!lead[j]) continue;
        sample[j] = (cur[k]^y[k]);
        ++k;

        // add sample[j] *column[j] to answer
        for(int i = 0; i < m; ++i) cur[i] ^= (sample[j]&M[i][j]);
    }

    // check answer wrt rref
    for(int i = 0; i < m; ++i) assert(cur[i] == y[i]); 


    return sample;
}

Vector sample_from_null(const Matrix &A, std::function<bit()> &coin) {

    Matrix M = rref(A);

    int m = M.size();
    assert(m); 
    int n = M[0].size();

    Vector y(m,0);

    Vector sample = get_random_solution(M,y,coin);
   // if(check_zero(sample)) cerr << "Warning: zero sample\n";
    return sample;
}

// only works if a solution exists otherwise assertion will fail
Vector sample_from_solutions(const Matrix &A, const Vector &y, std::function<bit()> &coin){
    Matrix M = A;

    int m = M.size();
    assert(m); 
    int n = M[0].size();
    
    // Append y as a column to M to perform rref on both
    for(int i = 0; i < m; ++i) M[i].push_back(y[i]);
    M = rref(M);

    // Remove the last column of M, then we need to solve for Mx = y1
    Vector y1(m);
    for(int i = 0; i < m; ++i ){
        y1[i] = M[i].back();
        M[i].pop_back();
    }

    Vector sample = get_random_solution(M,y1,coin);
    // check answer wrt original matrix
    Vector check = mul(A,sample);
    assert(check.size() == y.size());
    for(int i = 0; i < m; ++i) assert(check[i] == y[i]);

    return sample;
}

// returns all the solutions to Ax=0
Matrix get_entire_null(const Matrix &A){
    Matrix M = rref(A);

    int m = M.size();
    assert(m); 
    int n = M[0].size();

    vector<bool> lead = get_pivots(M);

    // count how many free variables
    int free = 0;
    for(int j = 0; j < n; ++j) free += (!lead[j]);
    
    // not too many free solutions
    assert(free < 30);
    Matrix solutions;

    solutions.resize(1<<free);
    // all possibilities for free variables
    for(int mask = 0; mask < (1<<free); ++mask){
        // fix the free variables
        Vector sample(n);
        Vector cur(m,0);
        for(int j = 0, k=0; j < n; ++j){
            if(lead[j]) {
                continue;
            }
            // value of jth variable is determined by mask
            sample[j] = !!(mask&(1<<k));
            ++k;
            // add sample[j] * column[j] to the answer
            for(int i = 0; i < m; ++i) cur[i] ^= (sample[j]&M[i][j]);
        }

        // calculate the other variables, we want cur to equal to y1
        for(int j =0,k=0 ;j < n; ++j){
            if(!lead[j]) continue;
            sample[j] = cur[k];
            ++k;
            // add sample[j] *column[j] to answer
            for(int i = 0; i < m; ++i) cur[i] ^= (sample[j]&M[i][j]);
        }

        // check answer wrt rref
        assert(check_zero(cur));
    
        // check answer wrt original matrix
        Vector check = mul(A,sample);
        assert(check_zero(check));

        solutions[mask] = sample;
    }

    return solutions;
}


Vector mul(const Matrix &M, const Vector &v){

    int i,j;
    Vector out(M.size(),0);

    assert(M[0].size() == v.size());
    for(j =0 ; j < v.size(); ++j){
        for(i = 0; i < M.size(); ++i) out[i] ^= (M[i][j]&v[j]); 
    }

    return out;
}

Matrix mul(const Matrix &A, const Matrix &B){
   int m1,n,m2;
   
   m1 = A.size();
   n = A[0].size();
   assert(B.size() == n);
   m2 = B[0].size();
   Matrix C(m1);
   for(int i = 0; i < m1; ++i) C[i].resize(m2);

   for(int i = 0; i < m1; ++i){
     for(int j = 0; j < m2; ++j){
        C[i][j] =0;
        for(int k = 0; k < n; ++k) C[i][j] ^= (A[i][k]&B[k][j]);
      }
    }
    
    return C;
}

Matrix transpose(const Matrix &A){
    int n,m;
    m = A.size();
    assert(m);
    n = A[0].size();
    Matrix T(n);

    for(int i = 0; i < n; ++i){
        T[i].resize(m);
        for(int j = 0; j < m; ++j){
            T[i][j] = A[j][i]; 
        }
    }

    return T;
}

Vector gen_random_vector(int n, std::function<bit()> &coin){
    Vector V(n);
    for(int i = 0; i < n; ++i){
        V[i] = coin();
    }
    return V;
}


bit dot(const Vector &a, const Vector &b) {
    assert(a.size() == b.size()); 
    bit prod = 0;
    for(int i = 0; i < a.size(); ++i) prod ^= (a[i]&b[i]);
    return prod;
}


Matrix gen_random_matrix(int m, int n, std::function<bit()> &coin){
    Matrix A(m);
    for(int i = 0; i < m; ++i){
        A[i].resize(n);
        for(int j = 0; j < n; ++j) {
            A[i][j] = coin();
        }
    }
    return A;
}

void print(const Matrix &M){
    int i,j;

    for(i = 0; i < M.size(); ++i){
        for(j = 0; j < M[i].size(); ++j){
            cout << M[i][j] << " ";
        }
        cout << "\n";
    }
}

void print(const Vector &v,bool orient){
    int i;

    for(int i = 0; i < v.size(); ++i) {
        cout << v[i] << (orient ? "\n" : " ");
    }
    cout << "\n";
}


