using namespace arma; using namespace std;
#include "armadillo"
#include <random>
#include <iostream>
#include <fstream>
#include <cstdlib>


int main()
{
    int N = 500; // number of agents
    double m0 = 2; // initial money
    double lambda = atof(argv[1]); // 0.25
    //double alpha = atof(argv[2]); 0.0
    //double gamma = atof(argv[3]); 0.5
    int n_sims = 1e3; // number of simulations
    int sims_done = 0;
    Transactions T = Transactions(transactions, m, m0, N, lambda)
    vec m_dist = T.make_m_array(N, m0);

    for(i = 0; i < 1e3; i++){
        cout << "Simulation " << sims_done << endl;
        int** transactions = T.make_trans_matrix(N);
        vec m = T.make_m_array(N, m0);
        vec m2 = T.do_trans(transactions, m, m0, N, lambda);

        // sum over results
        for (int i = 0; i < N; i++){
            m_dist(i) += m2(i);
        }

    sims_done += 1;
    }

    for (int i = 0; i < N; i++){
            m_dist(i) /= n_sims;
    }

    write_to_file(finalmoneydist,N);
}

