#include "armadillo"
#include "transactions.h"
#include <random>
#include <iostream>
#include <fstream>
using namespace arma; using namespace std;


int main()
{
    // SET INITIAL CONDITIONS
    int N = 1000; // number of agents
    double m0 = 2; // initial money
    double lambda = 0.25;//atof(argv[1]); // 0.25
    double alpha = 2.0;//atof(argv[2]); //0.0
    double gamma = 0.0; //atof(argv[3]); 0.5
    int n_sims = 1e3; // number of simulations
    int n_bins = int(m0*30/0.05);//income bins array

    // MAKE A HISTOGRAM VECTOR
    vec histbins = linspace(0,m0*30,n_bins+1);
    vec hist_total = zeros(n_bins);
    // make transactions happen between agents for the number of wanted experiments
    for(int j = 0; j < n_sims; j++){
        Transactions T(m0, N, lambda, gamma, alpha, n_bins,histbins);
        T.do_trans(1e7); // do_trans(transactions, m, N, lambda, gamma, alpha)
        uvec hist = T.getHistogram(histbins);
        cout << j+1 << endl;
        for(int i = 0; i < n_bins; i++){
            hist_total(i) += hist(i);
        }
    }
    // average over all all income distributions for all experiments
    Transactions T_hist(m0, N, lambda, gamma, alpha, n_bins, histbins);
    hist_total /= n_sims;
    // OUTPUT HISTOGRAM
    T_hist.write_to_file(hist_total);

}



