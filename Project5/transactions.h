#ifndef TRANSACTIONS_H
#define TRANSACTIONS_H
#include "armadillo"
#include <random>
#include <iostream>
#include <fstream>


using namespace arma;
using namespace std;

class Transactions{//(mat transactions, vec m,vec m0, int N, double lambda, double gamma,double alpha){ //hvorfor sendes det inn i klassen?
public:
    Transactions(double m0, int N, double lambda, double gamma, double alpha, int n_bins, vec hist_bins);
    void do_trans(int n_trans);
    void write_to_file(vec histogram);

    vec m;
    mat transactions_matrix;
    int m_N;
    double m_lambda;
    double m_alpha;
    double m_gamma;
    int bins;
    uvec getHistogram(vec linbins);
    vec bin_interval;
private:

};
#endif
