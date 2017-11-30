#include "transactions.h"

void Transactions::do_trans(int n_trans = 1e4){// burde ikke trenge m0 til aa sendes inn?

    // number of transactions: at least 10^7
    int trans_done = 0;

    // create a uniform random number distribution for creating epsilon
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> doubleRNG(0,1);
    double sum_ij;
    // perform 10^7 transactions ideally
    while(trans_done < n_trans){
        // choose random agents i and j
        int i = rand() % m_N; int j = rand() % m_N;
        if(i != j){
            double epsilon = doubleRNG(gen); //create a random double [0,1]
            sum_ij = m(i) + m(j);
            // money of the two agents are changed via a random transaction decided by epsilon
            m(i) = m_lambda*m(i) + epsilon*(1.0-m_lambda)*sum_ij;
            m(j) = m_lambda*m(j) + (1.0-epsilon)*(1.0-m_lambda)*sum_ij;
            trans_done += 1;
            // check if sum is equal
            if (sum_ij - m(i) - m(j) > 1e-8){
                cout << "transaction did not conserve money" << endl;
            }
        }
    }
}

uvec Transactions::getHistogram(vec linbins){
    //double maxmoney = m_N*m0;
    uvec histogram = hist(m, linbins);
    //cout << histogram << endl;
    return histogram;
}


void Transactions::write_to_file(vec histogram){
    ofstream myfile;
    myfile.open("hist.dat");
    for (int i = 0; i < bins; i++){
        myfile << histogram(i) << " " << bin_interval(i) << endl;
    }
    myfile.close();
}



Transactions::Transactions(double m0, int N, double lambda, double gamma, double alpha, int n_bins, vec hist_bins){
    m_N = N;
    m = zeros(m_N) + m0;
    m_alpha = alpha;
    m_gamma = gamma;
    m_lambda = lambda;
    transactions_matrix = zeros(N,N);
    bins = n_bins;
    bin_interval = hist_bins;
}
