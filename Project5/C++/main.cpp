#include "armadillo"
#include "transactions.h"
#include <random>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include "transactions.h"
#include "transactions.cpp"

using namespace arma; using namespace std;


int main(int argc, char* argv[])
{


    // SET INITIAL CONDITIONS
    int N = 1000; // number of agents
    double m0 = 100; // initial money
    double lambda = 0;//atof(argv[1]); // 0.25
    double alpha = 2.;//atof(argv[2]); //0.0
    double gamma = 4.; //atof(argv[3]); 0.5
    int n_sims = 1e3; // number of simulations
    int n_bins = int(5000);//income bins array
    int n_trans = 1e7;

    // MAKE A HISTOGRAM VECTOR
    vec histbins = linspace(0,5000,n_bins + 1);
    vec hist_total = zeros(n_bins);
    uvec hist;
    vec hist_TOT = zeros(n_bins);

    ofstream myfile;
    myfile.open("Ealpha2_lambda0_gamma_4_N1000.dat");

    int my_rank, numprocs, idum;

    Transactions T(m0, N, lambda, gamma, alpha, n_bins, histbins);

    for (int i = 0; i < 100; i++)
    {
         T.do_trans(n_trans);
    }

    cout << "equilib reached" << endl;

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    int no_intervalls = n_sims/numprocs;
    int myloop_begin = my_rank*no_intervalls + 1;
    int myloop_end = (my_rank+1)*no_intervalls;
    if ( (my_rank == numprocs-1) &&( myloop_end < n_sims) ) myloop_end = n_sims;
    MPI_Bcast (&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //own seed to individual processors
    idum = -1-my_rank;
    srand(idum);
    // random starting point

    // make transactions happen between agents for the number of wanted experiments
    for(int cycles = myloop_begin; cycles <= myloop_end; cycles++)
    {

        T.do_trans(n_trans); // do_trans(transactions, m, N, lambda, gamma, alpha)

        cout << cycles << endl;
        hist = T.getHistogram(histbins);

        for (int i = 0; i < n_bins; i++)
        {
            hist_TOT(i) += hist(i);
        }

    }
    for (int i = 0; i < n_bins; i++)
    {
       MPI_Reduce(&hist_TOT(i), &hist_total(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    if (my_rank == 0)
    {
        for (int i = 0; i < n_bins; i++)
        {
            myfile << hist_total(i)/double(n_sims) << " " << histbins(i) << endl;
        }
    }



    myfile.close();



    MPI_Finalize ();






}



