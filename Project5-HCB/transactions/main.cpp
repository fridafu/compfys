#include <iostream>
#include <fstream>
//#include <mpi.h>
#include "armadillo"

using namespace std;
using namespace arma;


int main(int argc, char *argv[])
{
    // agents N
    int N = 500;
    // start capital
    double m0 = 1000;
    vec agents = zeros(N)+m0;
    cout << agents << endl;
    int MCsteps = 100;
    int N_transactions = 1000;
    double p;
    // do MCcycles
    for (int i = 0; i<MCsteps; i++){
        // make transactions happen 
        for (int j=0; j<N_transactions; j++){
            int i_index = distribution(gen);
            int j_index = distribution(gen);
            double eps_ = eps(gen);
            // finding probability of interaction to take place
            if (agents(i_index) - agents(j_index) == 0){
                p = 1.
            }
            else{
                // put in expression for 
            }
            if (eps(gen)<p && (i_index != j_index)){
                // make a transaction happen
            }
        }
    }




    return 0;
}
