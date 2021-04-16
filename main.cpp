#include <algorithm>
#include <iostream>
#include <cassert>
#include <functional>
#include <random>
#include <chrono>
#include "matrix.h"

using namespace std;

void test_sample_from_null(int n, int m, function<bit()> &coin) {
    Matrix A = gen_random_matrix(m,n,coin);
    //  print(A);
    //  printf("\n rref \n");
    Matrix M = rref(A);
    //  print(M);
    Vector nul =  sample_from_null(A, coin);

    print(nul);
    Vector ans = mul(A,nul); 
    assert(ans.size() == m);
    assert(check_zero(ans));
}

// check that columns of matrix are orthogonal to each other
void test_orthogonal(Matrix &M){
    Matrix T = transpose(M);
    for(int i = 0; i < T.size(); ++i){
        for(int j = i; j < T.size(); ++j){
            assert(!dot(T[i],T[j]));
        }
    }
}

// s is taken to be 10^{n-1}, P_s is such that g=0
// start with all 1's column and then repeateadly add columns 
// that are orthogonal to previous ones (in the row null space)
Matrix create_P_s(int m, int n, function<bit()> &coin){
    Matrix M(m);
    for(int i = 0; i < m; ++i) M[i].push_back(1);

    //transpose
    Matrix T = transpose(M);

    for(int i = 1; i < n; ++i) {
        Vector col = sample_from_null(T,coin);

        assert(col.size() == m); 

        // add the new column
        T.emplace_back(col);
        for(int j = 0; j < m; ++j) M[j].push_back(col[j]); 
    }

    test_orthogonal(M);
    return M;
}

// generate a random P which contains a P_s with g = 0
Matrix create_P(int m, int n, function<bit()> &coin){
    int m1 = m/2;
    if(m1&1) --m1;

    Matrix P = create_P_s(m1,n,coin);

    for(int i = m1; i < m; ++i){
        Vector rv = gen_random_vector(n-1,coin);
        // rows is non orthogonal to s
        rv.insert(rv.begin(),0);

        P.push_back(rv);
    }
    return P;
}


// Generate samples of the form y(u,v) = \sum_{p\in P} p <u,p> <v,p>
Vector get_prover_sample(const  Matrix &P, function<bit()> &coin){
    int n,m;
    m = P.size();
    n = P[0].size();

    Vector u = gen_random_vector(n,coin);
    Vector v = gen_random_vector(n, coin);

    Vector y(n,0);

    for(int i = 0; i < m; ++i) {
        if(dot(P[i],u)&dot(P[i],v)) {
            cnot(P[i],y);
        }
    }

    return y;
}

string to_string(Vector &v){
    string s;
    for(int i = 0; i < v.size(); ++i) s.push_back('0'+v[i]);
    return s;
}

// How many are unique? How many are nonzero?
void check_prover_samples(int ns, const Matrix &P, function<bit()> &coin) {
    int i;

    // secret string
    Vector s(P[0].size(), 0);
    s[0]=1;

    vector<string> samples;
    int nz =0;
    // generate ns samples
    for(i = 0; i < ns; ++i){
        Vector y = get_prover_sample(P,coin);
        assert(!dot(y,s));
        nz += check_zero(y);
        samples.push_back(to_string(y));
    }

    sort(samples.begin(), samples.end());
    int uniq = unique(samples.begin(),samples.end())-samples.begin();

    cout << "total: " << ns << " zero: " << nz  << " unique: " << uniq << "\n";
}

void run_iqp_protocol(function<bit()> &coin){
    int n = 100, m = 50;
    Matrix P = create_P(m,n,coin);
    //print(P);
    check_prover_samples(2*m, P, coin);
}

int main(){

    std::default_random_engine generator;
    std::uniform_int_distribution<bit> distribution(0,1);
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    std::function<bit()> coin = std::bind ( distribution, generator );

    //test_sample_from_null(100,50,coin); 

    run_iqp_protocol(coin);


    return 0;
}
