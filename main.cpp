#include <algorithm>
#include <iostream>
#include <cassert>
#include <functional>
#include <random>
#include <chrono>
#include "matrix.h"

#define DETAILED 
#define OUTPUT_RANK

using namespace std;

void test_sample_from_null(int m, int n, function<bit()> &coin) {
    Matrix A = gen_random_matrix(m,n,coin);
    //  print(A);
    //  printf("\n rref \n");
    //  Matrix M = rref(A);
    //  print(M);
    Vector nul =  sample_from_null(A, coin);

    print(nul);
    Vector ans = mul(A,nul); 
    assert(ans.size() == m);
    assert(check_zero(ans));
}

void test_sample_from_solutions(int m, int n, function<bit()> &coin) {
    Matrix A = gen_random_matrix(m,n,coin);
 //   print(A);

    // To ensure system has a solution
    Vector x = gen_random_vector(n,coin);
    Vector y = mul(A,x);

//   print(y);
    Vector sol =  sample_from_solutions(A, y, coin);
//    print(sol);
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

// check that M^TM is symmetric and rank 1
void test_g1(Matrix &M){
    int i,j,n;
    Matrix A =mul(transpose(M),M);
    assert(rnk(A) == 1);
    n = A.size();
    assert(A[0].size() == n);

    for(i = 0; i < n; ++i) {
      for(j = 0; j < n; ++j){
        assert(A[i][j] == A[j][i]);  
      }
    }
}

// s is taken to be 10^{n-1}, Creates P_s such that g=1
// start with all 1's column
// let A = P_s^TP_s, then A[i][0] = A[i][i] 
// If A[i][i] = 0, then A[j][i] = A[i][j] = 0 for all j
Matrix create_P_s1(int m, int n, function<bit()> &coin) {
    Matrix M(m);
    for(int i = 0; i < m; ++i) M[i].push_back(1);
    // the nonzero row of A;
    Vector nz_row; 
    nz_row.push_back(1);

    Matrix T = transpose(M);
    
    for(int i = 1; i< n; ++i) {
       Vector col;
       
       // The ith row of A (A=M^TM), is either all zeros or nz_row, pick randomly
       if(coin()) {
          col = sample_from_solutions(T,nz_row, coin);
          nz_row.push_back(1);
       } else{
          col = sample_from_null(T,coin);
          nz_row.push_back(0);
       }

       assert(col.size() == m); 
       T.emplace_back(col);
       for(int j = 0; j < m; ++j) M[j].push_back(col[j]); 
    }
    
    test_g1(M);
#ifdef OUTPUT_RANK
    cout<< rnk(M) << ", ";
#endif
    return M;
}

// s is taken to be 10^{n-1}, Creates P_s such that g=0
// start with all 1's column and then repeateadly add columns 
// that are orthogonal to previous ones (in the row null space)
Matrix create_P_s0(int m, int n, function<bit()> &coin){
    Matrix M(m);
    for(int i = 0; i < m; ++i) M[i].push_back(1);

    //transpose
    Matrix T = transpose(M);

    int zero=0;
    for(int i = 1; i < n; ++i) {
        Vector col = sample_from_null(T,coin);

        assert(col.size() == m); 

        // add the new column
        T.emplace_back(col);
        for(int j = 0; j < m; ++j) M[j].push_back(col[j]); 
        if(check_zero(col)) ++zero;
    }

#ifdef OUTPUT_RANK
    cout<< rnk(M) << ", ";
#endif
    //cout << "zero " << zero << " rank of P_s " << rnk(M) << "\n";
    test_orthogonal(M);
    return M;
}

// generate a random P which contains a P_s with g as specified (0 or 1)
Matrix create_P(int g, int m, int n, function<bit()> &coin){
    int m1 = m/2;

    Matrix P;
    if(g == 0){
        // m1 has to be even for g=0
        if(m1&1) --m1;
        P = create_P_s0(m1,n,coin);
    } else {
        if(!(m1&1)) --m1;
        P = create_P_s1(m1,n,coin);
    }

    Matrix R;
    for(int i = m1; i < m; ++i){
        Vector rv = gen_random_vector(n-1,coin);
        // rows is non orthogonal to s
        rv.insert(rv.begin(),0);

        P.push_back(rv);
        R.push_back(rv);
    }

#ifdef OUTPUT_RANK
    cout << rnk(R) <<  ", " << rnk(P) << ", ";
#endif
    return P;
}

// Generate samples of the form y(u,v) = \sum_{p\in P} p <u,p> <v,p>
Vector get_prover_sample(const  Matrix &P, Vector &u, Vector &v){
    int n,m;
    m = P.size();
    n = P[0].size();

    Vector y(n,0);

    for(int i = 0; i < m; ++i) {
        if(dot(P[i],u)&dot(P[i],v)) {
            cnot(P[i],y);
        }
    }

    return y;
}


Matrix extract_P_s(const Matrix &P, const Vector &s){
    int n,m;
    m = P.size();

    Matrix P_s;

    for(int i = 0; i < m; ++i) {
        if(dot(P[i],s))  P_s.push_back(P[i]);
    }

    return P_s;
}

string to_string(Vector &v){
    string s;
    for(int i = 0; i < v.size(); ++i) s.push_back('0'+v[i]);
    return s;
}

// strategy for g = 0
Matrix get_samplesg0(int ns, const Matrix &P, function<bit()> &coin) {
    int n = P[0].size(); 
    Matrix Q(ns);

    for(int i = 0; i < ns; ++i){
        Vector u = gen_random_vector(n,coin);
        Vector v = gen_random_vector(n, coin);
        Vector y = get_prover_sample(P,u,v);
        Q[i] = y;
    }
   
   return Q;
}


Matrix get_samplesg1(int ns, const Matrix &P, function<bit()> &coin) {
    int n = P[0].size();
    Matrix Q(ns);
    Vector u = gen_random_vector(n,coin);
    // strategy for g = 1 uses the fixed u and random v each time
    for(int i = 0; i < ns; ++i){
        Vector v = gen_random_vector(n, coin);
        Vector y = get_prover_sample(P,u,v);
        Q[i] = y;
    }

    return Q;
}


Matrix get_samplesg1_basis(const Matrix &P, function<bit()> &coin) {
    int n = P[0].size(); 
    Matrix Q(n);
    Vector u = gen_random_vector(n,coin);
    
    for(int i = 0; i < n; ++i) {
        Vector v(n,0);
        v[i] = 1;
        Vector y = get_prover_sample(P,u,v);
        Q[i] = y;
    }

    return Q;
}

// How many are unique? How many are nonzero? How many are orthogonal to s
void check_prover_samples(int g, int ns, const Matrix &P, function<bit()> &coin) {
    // secret string
    int n = P[0].size();
    Vector s(n, 0);
    s[0]=1;

    vector<string> samples;
    int nz =0;
    // generate ns samples
    Matrix Q;
    //get appropriate samples;
    if(g == 0) Q = get_samplesg0(ns,P,coin);
    if(g == 1) Q = get_samplesg1(ns,P,coin);

    int orth=0;
    for(int i = 0; i < ns; ++i){
        orth += !dot(Q[i],s);
        nz += check_zero(Q[i]);
        samples.push_back(to_string(Q[i]));
    }

    // for g = 0 all samples are orthogonal
    if(g == 0) assert(orth == ns);

    sort(samples.begin(), samples.end());
    int uniq = unique(samples.begin(),samples.end())-samples.begin();
//    cout << "total: " << ns << " orth: " << orth << " zero: " << nz  << " unique: " << uniq-nz << "\n";
#ifdef OUTPUT_RANK
    cout << rnk(Q) << "\n";
#endif
}


long long estimate_search_space(const Matrix &P, function<bit()> &coin){
    int m = P.size();
    int n = P[0].size();
    
    Vector s(n, 0);
    s[0]=1;
    // total number of strings tried
    long long tries =0;   
    int attempts = 0;

    int rnkP = rnk(P);

    while(1) {
        attempts++;
        Matrix Q = get_samplesg1_basis(P, coin);
        int r = rnk(Q);

        tries += 1ll<<(rnkP-r);

        int all_orth = 1;
        for(int i = 0; i < Q.size(); ++i){
            all_orth &= !dot(Q[i],s);
        }

        // s lies in null of Q
        if(all_orth) break;
    }
    
    cout << attempts << ", ";
    return tries;
}

Matrix reduce_column(const Matrix &P){
   Matrix T = rref(transpose(P));

   while(!T.empty() && check_zero(T.back())){
    T.pop_back();
   }
   
   #ifdef OUTPUT_RANK
   cout << rnk(T) << ", " << T.size() << "\n";
   #endif
   return transpose(T);
}

// check that scrambled s is correct for a scrambled P
bool check_secret_1(const Matrix &P, const Vector &s){
    // check that first m_s rows are orthogonal to s
    int m = P.size(); 
    int m_s = m/2;
    if(!(m_s&1)) --m_s;
    
    // first m_s rows are non orthogonal to s
    bool ret = true;
    int c1=0,c2=0;
    for(int i = 0; i < m_s; ++i){
        if(dot(P[i],s) == 0) ret = false, ++c1; 
    }

    // last m-m_s rows are orthogonal to s
    for(int i = m_s; i < m; ++i){
        if(dot(P[i],s) == 1) ret = false, ++c2;
    } 

    if(!ret) {
        cout << "wrong at " << c1 << "/" << m_s << " "  << c2 << "/" << m-m_s << "  positions\n";
    
        for(int i = 0; i < m; ++i) {
            if(dot(P[i],s)==1) print(P[i]);
        }
    }

    printf("\n");

    return ret;
}

Matrix generate_Q1(const Matrix &P, function<bit()> &coin, int maxk){
    Matrix Q = get_samplesg1_basis(P, coin);
        
    for(int k = 1; k < maxk; ++k){
        Matrix Q2 = get_samplesg1_basis(P, coin);
        for(int i =0 ; i < Q2.size(); ++i) Q.emplace_back(Q2[i]);
    }
    
    return Q;
}

int find_secret_1(const Matrix &P, function<bit()> &coin){
    int m = P.size();
    int n = P[0].size();
    int attempts = max(m*m,500);
    
    if(rnk(P) != n) {
        cout << "-1, ";
        return 0;
    }

    cout << "RATIO " << " " << n*1.0/m << "\n";
    
    Vector s(n, 0);
    s[0]=1;
    // total number of strings tried
    int tries =0;
    while(attempts--) {
        Matrix Q = get_samplesg1_basis(P, coin);
        
        for(int k = 0; (1<<(1<<k)) < m; ++k){
            Matrix Q2 = get_samplesg1_basis(P, coin);
            for(int i =0 ; i < Q2.size(); ++i) Q.emplace_back(Q2[i]);
        }

        int r = rnk(Q);
        // assume all of Q is orthogonal to S.
        #ifdef DETAILED
        cout << "rank of Q: " << r << "\n";
        #endif
        if(n-r > 20) {
            #ifdef DETAILED
            cout << "Rank too low to find s\n";
            #endif
            continue;
        }
        
        Matrix trial = get_entire_null(Q);
        
        #ifdef DETAILED
        cout << trial.size() << " candidates to try\n";
        #endif
        random_shuffle(trial.begin(),trial.end());
        for(int i = 0; i < trial.size(); ++i){
            ++tries;
            Matrix P_s = extract_P_s(P,trial[i]);
                        
            // special case
            if(P_s.empty()) continue;

            Matrix A =mul(transpose(P_s),P_s);
            if(rnk(A) == 1) {
                #ifdef DETAILED
                cout << "Found\n";
                #endif
               
                if(!check_secret_1(P,trial[i])) {
                   #ifdef DETAILED
                   cout << "Warning: wrong string\n";
                   #endif
                   cout << "2, ";
                   continue;
                } else {
                   cout << 0 <<", ";
                }

                return tries;
            }
        }

    }

   #ifdef DETAILED
   cerr << "Reached max attempts\n";
   #endif
   cout << "1, ";
   return tries;
}


void run_iqp_protocol(int g, function<bit()> &coin){

    int n = 50, m = 50;
    cin >> m >> n;
    Matrix P = create_P(g,m,n,coin);

    vector<string> S;
    for(int i = 0; i < m; ++i) S.push_back(to_string(P[i]));
    for(int i = 0; i < m; ++i) for(int j = i+1; j < m; ++j) {
        if(S[i] == S[j]) {
            cout << S[i] << "\n" << S[j] << "\n\n";
        }
    }
    
    Matrix R = reduce_column(P);
   // int r = rnk(P);
   // for(int k = 1; k <= 10; ++k){

   //     Matrix Q = generate_Q1(R,coin,k);
   //     cout << r-rnk(Q) << ", ";
   // }
   // cout << r << "\n";
    // check_prover_samples(g,m*m, P, coin);
    // if(g == 1){ 
    //   auto tries = find_secret_1(R,coin);
    //    cout << tries <<"\n";
    //}
}

// rank of P^TP for a random P
void random_P_rank(function<bit()> &coin) {
   int n = 50, m = 50;
   cin >>n;
   m = n;
   Matrix M = gen_random_matrix(m,n,coin);
   Matrix A = mul(transpose(M),M);
   cout << rnk(A) << "\n";
}

int main(){

    std::default_random_engine generator;
    std::uniform_int_distribution<bit> distribution(0,1);
    auto seed= std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(seed);
    std::function<bit()> coin = std::bind ( distribution, generator );
    cout << seed << ", ";
    //test_sample_from_null(100,50,coin); 
   // test_sample_from_solutions(300,1000,coin); 

    run_iqp_protocol(1,coin);

    //random_P_rank(coin);
    return 0;
}
