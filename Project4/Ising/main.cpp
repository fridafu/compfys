#include <iostream>
#include "armadillo"
#include "ising.h"
#include <iostream>
#include <fstream>
#include <mpi.h>

using namespace std;
using namespace arma;


int main(int nargs, char* args[])
{

    double k = 1.38064852e-23;
    double J = 1;
    double That = 1;
    double T = abs(That*J/k);


    //b)
    /*
    //
    int L = 2;
    double beta = 1./(T*k);
    int steps = 1000000;
    double z = 2*(exp(-8*J*beta) + exp(8*J*beta)) + 12;
    double M  = 0;
    double JB8 = 8*beta*J;
    double mE = -8*J*(sinh(JB8)/(cosh(JB8)+3));
    double mE2 = 64*J*J*(cosh(JB8)/(cosh(JB8)+3));
    double mM2 = 8*((exp(JB8) + 1)/(cosh(JB8) + 3)) ;
    double absM = 2*((exp(JB8) + 2)/(cosh(JB8) + 3)) ;
    double Cv = ((-128*J*J)/(T*T))*((1/z)*cosh(JB8) - (8/(z*z))*sinh(JB8)*sinh(JB8));
    double X = (8*(exp(JB8) + 1)/(cosh(JB8) + 3) - (2*(exp(JB8) + 2)/(cosh(JB8) + 3))*(2*(exp(JB8) + 2)/(cosh(JB8) + 3)))/(k*T);

    Ising L1 = Ising(J,L,T);
    L1.exp_vals(steps);
    vec expvals = L1.get_expectation_values();
    cout << "calculated energy = " << expvals(0) << endl;
    cout << "analytical energy = " << mE << endl;
    cout << "calculated energy^2 = " << expvals(1) << endl;
    cout << "analytical energy^2 = " << mE2 << endl;
    cout << "calculated magnetization = " << expvals(2) << endl;
    cout << "analytical magnetization = " << M << endl;
    cout << "calculated abs(magnetization) = " << expvals(4) << endl;
    cout << "analytical abs(magnetization) = " << absM << endl;
    cout << "calculated magnetization^2 = " << expvals(3) << endl;
    cout << "analytical magnetization^2 = " << mM2 << endl;
    cout << "calculated Cv = " << L1.heat_capacity() << endl;
    cout << "analytical Cv = " << Cv << endl;
    cout << "calculated X = " << L1.magnetic_susceptibility() << endl;
    cout << "analytical X = " << X << endl;

    cout << "_______________________________________________" << endl;
    L1.reset_expectation_values();
    L1.exp_vals(10*steps);

    vec expvals2 = L1.get_expectation_values();
    cout << "HERE IS ANOTHER ROUND" << endl;
    cout << "calculated energy = " << expvals2(0) << endl;
    cout << "analytical energy = " << mE << endl;
    cout << "calculated energy^2 = " << expvals2(1) << endl;
    cout << "analytical energy^2 = " << mE2 << endl;
    cout << "calculated magnetization = " << expvals2(2) << endl;
    cout << "analytical magnetization = " << M << endl;
    cout << "calculated abs(magnetization) = " << expvals2(4) << endl;
    cout << "analytical abs(magnetization) = " << absM << endl;
    cout << "calculated magnetization^2 = " << expvals2(3) << endl;
    cout << "analytical magnetization^2 = " << mM2 << endl;
    cout << "calculated Cv = " << L1.heat_capacity() << endl;
    cout << "analytical Cv = " << Cv << endl;
    cout << "calculated X = " << L1.magnetic_susceptibility() << endl;
    cout << "analytical X = " << X << endl;


    /*
terminal output:
calculated energy = -7.98461
analytical energy = -7.98393
calculated energy^2 = 63.8769
analytical energy^2 = 63.8714
calculated magnetization = 0.077706
analytical magnetization = 0
calculated abs(magnetization) = 3.99494
analytical abs(magnetization) = 3.99464
calculated magnetization^2 = 15.9745
analytical magnetization^2 = 15.9732
calculated Cv = 1.6968e-24
analytical Cv = 6.06309e-45
calculated X = 0.0149544
analytical X = 0.016043
_______________________________________________
HERE IS ANOTHER ROUND
calculated energy = -7.9839
analytical energy = -7.98393
calculated energy^2 = 63.8712
analytical energy^2 = 63.8714
calculated magnetization = -0.0677004
analytical magnetization = 0
calculated abs(magnetization) = 3.99461
analytical abs(magnetization) = 3.99464
calculated magnetization^2 = 15.9731
analytical magnetization^2 = 15.9732
calculated Cv = 1.77514e-24
analytical Cv = 6.06309e-45
calculated X = 0.0162285
analytical X = 0.016043
    */

    //c)
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


    //accepted configurations
    /*
    int LL = 20;
    int acceptconfig1;
    int acceptconfig2;

    vec temperature = linspace(1*T, 2.4*T, 20);
    vec that = linspace(1,2.4,20);
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
        onesconfig.open("acceptconfigones" + to_string(j+1) + ".txt");
        randomconfig.open("acceptconfigrand" + to_string(j+1) + ".txt");
        for (int i = 0; i < totsteps/checkstep; i++)
        {



            L1.exp_vals(checkstep);
            acceptconfig1 = L1.get_configurations();
            L2.exp_vals(checkstep);
            acceptconfig2 = L2.get_configurations();

            onesconfig << acceptconfig1 << " " << that(j) << " " << (i+1)*checkstep << endl;
            randomconfig << acceptconfig2 << " " << that(j) << " " << (i+1)*checkstep << endl;


        }
        onesconfig.close();
        randomconfig.close();
    {
    */

    //d)
    /*
    int L = 20;

    int initsteps = 10000;
    int steps = 100000;
    vec exp_vals;
    double PE;
    vec energy = linspace(-820,-750, 71);

    for (int i = 0; i < size(energy)(0); i++)
    {
        Ising L1 = Ising(J,L,T);


        PE = L1.energy_probability(initsteps,steps,energy(i));
        exp_vals = L1.get_expectation_values();
        if (PE != 0)
        {
            cout << "P(" << energy(i) << ") = " << PE << endl;
            cout << "Variance sigma_squared = " << exp_vals(1) - exp_vals(0)*exp_vals(0) << endl;
            cout << "__________________________________________________________________" << endl;
         }
    }
    */


    /*
Terminal output:
P(-800) = 0.8698
Variance sigma_squared = 9.31287
__________________________________________________________________
P(-792) = 0.11579
Variance sigma_squared = 9.32761
__________________________________________________________________
P(-788) = 0.00451
Variance sigma_squared = 9.32324
__________________________________________________________________
P(-784) = 0.00781
Variance sigma_squared = 9.3839
__________________________________________________________________
P(-780) = 0.0004
Variance sigma_squared = 9.14445
__________________________________________________________________
P(-776) = 0.00046
Variance sigma_squared = 9.51847
__________________________________________________________________
P(-772) = 0.0001
Variance sigma_squared = 9.46302
__________________________________________________________________
P(-768) = 1e-05
Variance sigma_squared = 9.24947
__________________________________________________________________
P(-764) = 1e-05
Variance sigma_squared = 9.4708
__________________________________________________________________
*/
    char *outfilename;
    long idum;
    int **spin_matrix, n_spins, mcs, my_rank, numprocs; double w[17], average[5], total_average[5],
        initial_temp, final_temp, E, M, temp_step;


    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0 && argc <= 1)
    {
        cout << "Bad Usage: " << argv[0] << " read output file" << endl;
        exit(1);
    }

    n_spins = 40; mcs = 1000000; initial_temp = 2.4; final_temp = 2.7; temp_step =0.1;

    MPI_Bcast (&n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // every node has its own seed for the random numbers, this is important else
    // if one starts with the same seed, one ends with the same random numbers
    idum = -1-my_rank;
    // random starting point

    //make loop over temperatures

    //perform monte carlo here

    //end monte carlo here

    //find total average

    for( int i =0; i < 5; i++)
    {
        MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    //end loop over temperatures here

    MPI_Finalize ();



    return 0;

}
