#include <iostream>
#include "armadillo"
#include "ising.h"

using namespace std;
using namespace arma;


int main()
{
    double k = 1.38064852e-23;
    double J = -1;
    double That = 1;
    double T = That*J/k;
    int L = 2;
    double beta = 1/(k*T);
    int steps = 100000;


    Ising L2 = Ising(J,L,T);
    vec exp_v = L2.exp_vals(steps);
    cout << exp_v << endl;
    cout << "calculated CV = " << L2.heat_capacity() << endl;
    cout << "calculated X = " << L2.magnetic_susceptibility() << endl;
    double z = 12 + 2*(exp(-8*beta*J) + exp(8*beta*J));
    double mE = (16*J*(exp(-8*beta*J) + exp(8*beta*J)))/z;
    double mE2 = 128*J*J*(exp(-8*beta*J) + exp(8*beta*J))/z;
    double mM2 = (32*exp(8*beta*J) + 32)/z;
    double absM = (8/z)*(exp(8*beta*J) + 2);
    double cv = (J*J/(k*T*T))*(mE2 - mE*mE);
    double magsus = (1/(k*T))*(mM2 - absM*absM);

    cout << "Analytic E = " << mE << endl;
    cout << "Analytic E2 = " << mE2 << endl;

    cout << "analytical CV = " << cv << endl;
    cout << "analytical X = " << magsus << endl;

    return 0;
}
