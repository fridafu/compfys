#include <iostream>
#include "armadillo"
#include "ising.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;


int main()
{
    double k = 1.38064852e-23;
    double J = 1;
    double That = 1;
    double T = abs(That*J/k);
    //b)


    int L = 2;
    double beta = 1./(T*k);
    int steps = 10000000;
    double z = 2*(exp(-8*J*beta) + exp(8*J*beta)) + 12;
    double Cv = ((-128*J*J)/(T*T))*((1/z)*cosh(8*beta*J) - (8/(z*z))*sinh(8*beta*J)*sinh(8*beta*J));
    double X = (8*(exp(8*beta*J) + 1)/(cosh(8*beta*J) + 3) - (2*(exp(8*beta*J) + 2)/(cosh(8*beta*J) + 3))*(2*(exp(8*beta*J) + 2)/(cosh(8*beta*J) + 3)))/(k*T);

    Ising L1 = Ising(J,L,T);
    L1.exp_vals(steps);
    vec expvals = L1.get_expectation_values();

    cout << "calculated energy = " << expvals(0) << endl;
    cout << "calculated Cv = " << L1.heat_capacity() << endl;
    cout << "calculated X = " << L1.magnetic_susceptibility() << endl;
    cout << "analytical Cv = " << Cv << endl;
    cout << "analytical X = " << X << endl;



    //b)

    ofstream myfileones;
    ofstream myfilerandom;
    int LL = 20;
    Ising L2 = Ising(J,LL,T);
    mat initstate1 = ones(LL,LL);
    L2.set_state(initstate1);

    Ising L3 = Ising(J,LL,T);

    int totsteps = 5000;
    int checkstep = 10;
    vec expvals1;
    vec expvals2;

    myfileones.open("expvalsones.txt");
    myfilerandom.open("expvalsrandom.txt");
    for (int i = 0; i < totsteps/checkstep; i++)
    {

        L2.exp_vals(checkstep);
        expvals1 = L2.get_expectation_values();
        L3.exp_vals(checkstep);
        expvals2 = L3.get_expectation_values();
        myfileones << expvals1(0) << " " << expvals1(4) << " " << (i+1)*checkstep << endl;
        myfilerandom << expvals2(0) << " " << expvals2(4) << " " << (i+1)*checkstep << endl;
    }
    myfileones.close();
    myfilerandom.close();

    return 0;
}
