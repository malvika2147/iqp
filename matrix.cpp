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

Vector sample_from_null(const Matrix &A, std::function<bit()> &coin) {

    Matrix M = rref(A);

    int m = M.size();
    assert(m); 
    int n = M[0].size();
    vector<bool> lead(n, false);

    Vector sample(n);
    Vector cur(m,0);
    
    int rank=0;
    for(int i = 0,j=0; i < m; ++i,++j) {
        while(j < n && !M[i][j]) ++j;
        if(j >= n) break;

        lead[j] = 1;
        ++rank;
    } 

    // fix the free variables
    for(int j = 0; j < n; ++j){
        if(lead[j]) {
            continue;
        }
        sample[j] = coin();
        // add sample[j] * column[j] tp the answer
        for(int i = 0; i < m; ++i) cur[i] ^= (sample[j]&M[i][j]);
    }

    // calculate the other variables
    for(int j =0,k=0 ;j < n; ++j){
        if(!lead[j]) continue;
        sample[j] = cur[k++];

        // add sample[j] *column[j] to answer
        for(int i = 0; i < m; ++i) cur[i] ^= (sample[j]&M[i][j]);
    }

    for(int i = 0; i < m; ++i) assert(!cur[i]); 

    return sample;
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


