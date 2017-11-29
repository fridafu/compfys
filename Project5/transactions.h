#ifndef TRANSACTIONS_H
#define TRANSACTIONS_H
#include "armadillo"
#include <random>
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace arma;
using namespace std;

class Transactions(transactions, m, m0, N, lambda, gamma, alpha)
{
public:

    vec do_trans(transactions, m, N, lambda, gamma, alpha, m0);

    vec make_m_array(N, m_0);

    void make_trans_matrix(N);

    void write_to_file(N);

private:

}
