#include <iostream>
#include "armadillo"
#include "ising.h"

using namespace std;
using namespace arma;


int main()
{
    double k = 1.38064852e-23;
    double J = 1;
    double That = 1;

    int L = 2;
    double beta = 1./(That);
    int steps = 10000000;


    Ising L2 = Ising(J,L,That);
    L2.exp_vals(steps);

    vec expvals = L2.get_expectation_values();
    cout << "calculated energy = " << expvals(0) << endl;
    cout << " calculated Cv = " << L2.heat_capacity();



    return 0;
}
