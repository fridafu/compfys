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
    double T = abs(That*J/k);

    int L = 2;
    double beta = 1./(T*k);
    int steps = 10000000;


    Ising L2 = Ising(J,L,T);
    L2.exp_vals(steps);

    vec expvals = L2.get_expectation_values();
    cout << "calculated energy = " << expvals(0) << endl;
    cout << "calculated Cv = " << L2.heat_capacity() << endl;
    cout << "calculated X = " << L2.magnetic_susceptibility() << endl;


    double z = 2*(exp(-8*J*beta) + exp(8*J*beta)) + 12;
    double Cv = ((-128*J*J)/(T*T))*((1/z)*cosh(8*beta*J) - (8/(z*z))*sinh(8*beta*J)*sinh(8*beta*J));

    cout << "analytical Cv = " << Cv << endl;

    double X = (8*(exp(8*beta*J) + 1)/(cosh(8*beta*J) + 3) - (2*(exp(8*beta*J) + 2)/(cosh(8*beta*J) + 3))*(2*(exp(8*beta*J) + 2)/(cosh(8*beta*J) + 3)))/(k*T);

    cout << "analytical X = " << X << endl;



    return 0;
}
