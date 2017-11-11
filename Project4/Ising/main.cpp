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

    /*
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

    */

    //b)
    /*
    int LL = 20;
    ofstream myfileones;
    ofstream myfilerandom;
    ofstream myfileones24;
    ofstream myfilerandom24;
    mat initstate1 = ones(LL,LL);

    Ising L1 = Ising(J,LL,T);
    L1.set_state(initstate1);
    Ising L2 = Ising(J,LL,T);
    Ising L3 = Ising(J,LL,2.4*T);
    L3.set_state(initstate1);
    Ising L4 = Ising(J,LL,2.4*T);





    int totsteps = 10000;
    int checkstep = 10;
    vec expvals1;
    vec expvals2;
    vec expvals3;
    vec expvals4;

    myfileones.open("expvalsones.txt");
    myfilerandom.open("expvalsrandom.txt");
    myfileones24.open("expvalsones24.txt");
    myfilerandom24.open("expvalsrandom24.txt");
    for (int i = 0; i < totsteps/checkstep; i++)
    {

        L1.exp_vals(checkstep);
        expvals1 = L1.get_expectation_values();
        L2.exp_vals(checkstep);
        expvals2 = L2.get_expectation_values();
        L3.exp_vals(checkstep);
        expvals3 = L3.get_expectation_values();
        L4.exp_vals(checkstep);
        expvals4 = L4.get_expectation_values();

        myfileones << expvals1(0) << " " << expvals1(4) << " " << (i+1)*checkstep << endl;
        myfilerandom << expvals2(0) << " " << expvals2(4) << " " << (i+1)*checkstep << endl;
        myfileones24 << expvals3(0) << " " << expvals3(4) << " " << (i+1)*checkstep << endl;
        myfilerandom24 << expvals4(0) << " " << expvals4(4) << " " << (i+1)*checkstep << endl;
    }
    myfileones.close();
    myfilerandom.close();
    myfileones24.close();
    myfilerandom24.close();
    */

    int LL = 20;
    int acceptconfig1;
    int acceptconfig2;

    vec temperature = linspace(1*T, 2.4*T, 10);
    mat initstate = ones(LL,LL);
    int totsteps = 10000;
    int checkstep = 10;




    for (int j = 0; j < size(temperature)(0); j++)
    {
        Ising L1 = Ising(J, LL, temperature(j));
        Ising L2 = Ising(J, LL, temperature(j));
        L1.set_state(initstate);
        ofstream onesconfig;
        ofstream randomconfig;
        for (int i = 0; i < totsteps/checkstep; i++)
        {

            onesconfig.open("acceptconfigones" + to_string(j+1) + ".txt");
            randomconfig.open("acceptconfigrand" + to_string(j+1) + ".txt");

            L1.exp_vals(checkstep);
            acceptconfig1 = L1.get_configurations();
            L2.exp_vals(checkstep);
            acceptconfig2 = L2.get_configurations();

            onesconfig << acceptconfig1 << " " << temperature(j) << " " << (i+1)*checkstep << endl;
            randomconfig << acceptconfig2 << " " << temperature(j) << " " << (i+1)*checkstep << endl;


        }
        onesconfig.close();
        randomconfig.close();

    }


    return 0;

}
