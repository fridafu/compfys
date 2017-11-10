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
    int steps = 10000;

    Ising L2 = Ising(J,L,T);
    L2.exp_vals(steps);
    vec expvals = L2.get_expectation_values();
    cout << "calculated energy = " << expvals(0) << endl;
    cout << "calculated CV = " << L2.heat_capacity() << endl;
    cout << "calculated X = " << L2.magnetic_susceptibility() << endl;


    return 0;
}
