

void do_trans(transactions, m, N, lambda, gamma, alpha, m0){
    int n_trans = 1e4; 	// number of transactions: at least 10^7
    int trans_done = 0; int sims_done = 0;

    // create a uniform random number distribution for creating epsilon
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> doubleRNG(0,1);

    // perform 10^3 simulations
    while(trans_done < n_trans){
        // choose random agents i and j
        i = rand() % N; j = rand() % N;
        if(i != j){
            double epsilon = doubleRNG(gen); //create a random double [0,1]
            sum_ij = m(i) + m(j);
            // money of the two agents are changed via a random transaction decided by epsilon
            m(i) = lambda*m(i) + epsilon*(1.0-lambda)*sum_ij;
            m(j) = lambda*m(j) + (1.0-epsilon)*(1.0-lambda)*sum_ij;
            sims_done += 1;
        }
    }
    return m;
}


void write_to_file(vec m, int N){
    ofstream myfile;
    myfile.open("m.dat");
    for (int i = 0; i < N; i++){
        myfile << m(i) << endl;
    }
    myfile.close();
}

vec make_m_array(int N, int m0){
    vec m = new vec(N);
    for (int i = 0; i < N; i++){
        m(i) = m0;
    }
    return m;
}

int** make_trans_matrix(int N){
    int** transactions = new int*[N];
    for (int i = 0; i < N; i++){
        transactions[i] = new int[N];
    }
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            transactions[i][j] = 0;
        }
    }
    return transactions;
}
