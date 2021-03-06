#include "transactions.h"



void Transactions::do_trans(int n_trans = 1e4){// burde ikke trenge m0 til aa sendes inn?

    // number of transactions: at least 10^7
    int trans_done = 0;

    // create a uniform random number distribution for creating epsilon
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> doubleRNG(0,1);
    double sum_ij;
    double p_ij;
    double norm = 1/pow(400, m_alpha);
    int max_cij = 1;
    double max_pij = 1;
    double c_ij;
    double epsilon;
    int i;
    int j;
    double r;
    double dm = 1./m_0;
    // perform 10^7 transactions ideally
    while(trans_done < n_trans)
    {
        // choose random agents i and j
        i = rand() % m_N; j = rand() % m_N;
        r = doubleRNG(gen);
        // add probability to do transactions with agants with similar funds

        if (m(i) != m(j))
        {
            p_ij =  1000*pow(fabs((m(i)-m(j))/dm),-m_alpha)*pow((c_ij +1)/double(max_cij),m_gamma);
        }
        else
        {
            p_ij = 1;
        }

        if(i != j && p_ij>r){

            transactions_matrix(i,j) = c_ij + 1; // update number of transactions done for the agent i & j
            transactions_matrix(j,i) = c_ij + 1;
            if (c_ij + 1 > max_cij)
            {
                max_cij = c_ij + 1;
            }

            epsilon = doubleRNG(gen); //create a random double [0,1]
            sum_ij = m(i) + m(j);
            // money of the two agents are changed via a random transaction decided by epsilon
            m(i) = m_lambda*m(i) + epsilon*(1.0-m_lambda)*sum_ij;
            m(j) = m_lambda*m(j) + (1.0-epsilon)*(1.0-m_lambda)*sum_ij;
            trans_done += 1;
        }
    }
}

uvec Transactions::getHistogram(vec linbins){
    //double maxmoney = m_N*m0;
    uvec histogram = hist(m, linbins);

    return histogram;
}


void Transactions::write_to_file(vec histogram){
    // write to file, including parameters for the run
    ofstream myfile;
    myfile.open("Dhist_alpha15_N500_lam0.dat");
    myfile << "#alpha=" << m_alpha << endl;
    myfile << "#lambda=" << m_lambda << endl;
    myfile << "#gamma=" << m_gamma << endl;

    for (int i = 0; i < bins; i++)
    {
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
    m_0 = m0;
}
