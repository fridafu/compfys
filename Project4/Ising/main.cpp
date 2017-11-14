#include <iostream>
#include "armadillo"
#include "ising.h"
#include <iostream>
#include <fstream>
#include <mpi.h>
#include "ising.cpp"

using namespace std;
using namespace arma;


int main(int argc, char* argv[])
{

    /*
    double k = 1.38064852e-23;
    double J = 1;
    double That = 1;
    double T = abs(That*J/k);
    */
    /*
    //b)
    int L = 2;

    ofstream myfile;
    myfile.open("compare.txt");
    int L = 2;
    Ising L1 = Ising(J,L,T);



    double beta = 1./(T*k);
    int steps = 1000;
    double z = 2*(exp(-8*J*beta) + exp(8*J*beta)) + 12;
    double M  = 0;
    double JB8 = 8*beta*J;
    double mE = -8*J*(sinh(JB8)/(cosh(JB8)+3));
    double mE2 = 64*J*J*(cosh(JB8)/(cosh(JB8)+3));
    double mM2 = 8*((exp(JB8) + 1)/(cosh(JB8) + 3)) ;
    double absM = 2*((exp(JB8) + 2)/(cosh(JB8) + 3)) ;
    double Cv = ((-128*J*J)/(T*T))*((1/z)*cosh(JB8) - (8/(z*z))*sinh(JB8)*sinh(JB8));
    double X = (8*(exp(JB8) + 1)/(cosh(JB8) + 3) - (2*(exp(JB8) + 2)/(cosh(JB8) + 3))*(2*(exp(JB8) + 2)/(cosh(JB8) + 3)))/(k*T);


    L1.exp_vals(steps);
    vec expvals = L1.get_expectation_values();



    cout << "1000000 steps:" << endl;
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
    cout << "reset expectation values. 10000000 new steps:" << endl;
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

    */
    /*
terminal output:
1000000 steps:
calculated energy = -7.98421
analytical energy = -7.98393
calculated energy^2 = 63.8737
analytical energy^2 = 63.8714
calculated magnetization = 0.0602226
analytical magnetization = 0
calculated abs(magnetization) = 3.99474
analytical abs(magnetization) = 3.99464
calculated magnetization^2 = 15.9737
analytical magnetization^2 = 15.9732
calculated Cv = 1.74081e-24
analytical Cv = 6.06309e-45
calculated X = 0.0157248
analytical X = 0.016043
_______________________________________________
reset expectation values. 10000000 new steps:
HERE IS ANOTHER ROUND
calculated energy = -7.98396
analytical energy = -7.98393
calculated energy^2 = 63.8717
analytical energy^2 = 63.8714
calculated magnetization = -0.00419586
analytical magnetization = 0
calculated abs(magnetization) = 3.99465
analytical abs(magnetization) = 3.99464
calculated magnetization^2 = 15.9733
analytical magnetization^2 = 15.9732
calculated Cv = 1.76816e-24
analytical Cv = 6.06309e-45
calculated X = 0.0160374
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

    int initsteps = 100000;
    int steps = 100000;
    vec exp_vals;
    double PE;
    vec energy = linspace(-900,-760, 36);
    double meanE = 0;
    double meanE2 = 0;
    double totP = 0;
    ofstream myfile;
    myfile.open("prob.txt");

    for (int i = 0; i < size(energy)(0); i++)
    {
        Ising L1 = Ising(J,L,T);


        PE = L1.energy_probability(initsteps,steps,energy(i));
        exp_vals = L1.get_expectation_values();
        cout << energy(i) << endl;
        if (PE != 0)
        {
            meanE += energy(i)*PE;
            meanE2 += energy(i)*energy(i)*PE;
            totP += PE;

            cout << "P(" << energy(i) << ") = " << PE << endl;
            cout << "Variance sigma_squared = " << exp_vals(1) - exp_vals(0)*exp_vals(0) << endl;
            cout << "__________________________________________________________________" << endl;
            myfile << energy(i) << " " << PE << endl;
         }
    }
    myfile.close();
    cout << meanE2 << " " << meanE << " " << totP << endl;
    cout << "Computed variance = " << meanE2 - meanE*meanE << endl;
    */


/*
 * T = 1
Terminal output:
P(-800) = 0.869417
Variance sigma_squared = 9.38907
__________________________________________________________________
P(-792) = 0.117183
Variance sigma_squared = 9.38939
__________________________________________________________________
P(-788) = 0.0043593
Variance sigma_squared = 9.38291
__________________________________________________________________
P(-784) = 0.008031
Variance sigma_squared = 9.39139
__________________________________________________________________
P(-780) = 0.0005821
Variance sigma_squared = 9.34394
__________________________________________________________________
P(-776) = 0.0003893
Variance sigma_squared = 9.28129
__________________________________________________________________
P(-772) = 3.57e-05
Variance sigma_squared = 9.27876
__________________________________________________________________
P(-768) = 1.51e-05
Variance sigma_squared = 9.27582
__________________________________________________________________
P(-764) = 3.3e-06
Variance sigma_squared = 9.34417
__________________________________________________________________
P(-760) = 1e-07
Variance sigma_squared = 9.27119
__________________________________________________________________
P(-756) = 2e-07
Variance sigma_squared = 9.28069
__________________________________________________________________
638196 -798.872 1.00002
Computed variance = -1.01526

_____________

T = 2.4
Terminal output:
P(-304) = 2e-05
Variance sigma_squared = 3299.64
__________________________________________________________________
P(-308) = 2e-05
Variance sigma_squared = 3275.05
__________________________________________________________________
P(-312) = 1e-05
Variance sigma_squared = 3253.82
__________________________________________________________________
P(-316) = 3e-05
Variance sigma_squared = 3207.91
__________________________________________________________________
P(-320) = 4e-05
Variance sigma_squared = 3262.02
__________________________________________________________________
P(-324) = 0.00014
Variance sigma_squared = 3184.88
__________________________________________________________________
P(-328) = 0.00011
Variance sigma_squared = 3283.19
__________________________________________________________________
P(-332) = 0.00017
Variance sigma_squared = 3167.67
__________________________________________________________________
P(-336) = 0.00036
Variance sigma_squared = 3184.22
__________________________________________________________________
P(-340) = 0.00021
Variance sigma_squared = 3235.51
__________________________________________________________________
P(-344) = 0.00041
Variance sigma_squared = 3291.87
__________________________________________________________________
P(-348) = 0.00041
Variance sigma_squared = 3201.48
__________________________________________________________________
P(-352) = 0.00066
Variance sigma_squared = 3214.02
__________________________________________________________________
P(-356) = 0.00101
Variance sigma_squared = 3299.02
__________________________________________________________________
P(-360) = 0.00108
Variance sigma_squared = 3213.85
__________________________________________________________________
P(-364) = 0.00142
Variance sigma_squared = 3274.84
__________________________________________________________________
P(-368) = 0.00174
Variance sigma_squared = 3245.65
__________________________________________________________________
P(-372) = 0.00208
Variance sigma_squared = 3243.67
__________________________________________________________________
P(-376) = 0.00266
Variance sigma_squared = 3234.94
__________________________________________________________________
P(-380) = 0.00252
Variance sigma_squared = 3221.89
__________________________________________________________________
P(-384) = 0.00368
Variance sigma_squared = 3215.49
__________________________________________________________________
P(-388) = 0.00449
Variance sigma_squared = 3309.92
__________________________________________________________________
P(-392) = 0.00466
Variance sigma_squared = 3223.98
__________________________________________________________________
P(-396) = 0.00585
Variance sigma_squared = 3272.13
__________________________________________________________________
P(-400) = 0.00727
Variance sigma_squared = 3284.67
__________________________________________________________________
P(-404) = 0.00844
Variance sigma_squared = 3217.17
__________________________________________________________________
P(-408) = 0.00928
Variance sigma_squared = 3168.06
__________________________________________________________________
P(-412) = 0.01108
Variance sigma_squared = 3190.92
__________________________________________________________________
P(-416) = 0.01214
Variance sigma_squared = 3164.8
__________________________________________________________________
P(-420) = 0.01246
Variance sigma_squared = 3237.36
__________________________________________________________________
P(-424) = 0.01424
Variance sigma_squared = 3321.33
__________________________________________________________________
P(-428) = 0.01635
Variance sigma_squared = 3248.7
__________________________________________________________________
P(-432) = 0.017
Variance sigma_squared = 3235.66
__________________________________________________________________
P(-436) = 0.01947
Variance sigma_squared = 3258.43
__________________________________________________________________
P(-440) = 0.0199
Variance sigma_squared = 3222.91
__________________________________________________________________
P(-444) = 0.02133
Variance sigma_squared = 3302.75
__________________________________________________________________
P(-448) = 0.022
Variance sigma_squared = 3193.35
__________________________________________________________________
P(-452) = 0.02294
Variance sigma_squared = 3229.99
__________________________________________________________________
P(-456) = 0.02383
Variance sigma_squared = 3213.21
__________________________________________________________________
P(-460) = 0.02453
Variance sigma_squared = 3256.99
__________________________________________________________________
P(-464) = 0.02805
Variance sigma_squared = 3258.1
__________________________________________________________________
P(-468) = 0.02747
Variance sigma_squared = 3191.68
__________________________________________________________________
P(-472) = 0.02738
Variance sigma_squared = 3232.26
__________________________________________________________________
P(-476) = 0.02677
Variance sigma_squared = 3243.25
__________________________________________________________________
P(-480) = 0.02918
Variance sigma_squared = 3099.14
__________________________________________________________________
P(-484) = 0.02722
Variance sigma_squared = 3258.75
__________________________________________________________________
P(-488) = 0.02757
Variance sigma_squared = 3233.13
__________________________________________________________________
P(-492) = 0.02775
Variance sigma_squared = 3218.13
__________________________________________________________________
P(-496) = 0.02733
Variance sigma_squared = 3265.61
__________________________________________________________________
P(-500) = 0.02646
Variance sigma_squared = 3280.46
__________________________________________________________________
P(-504) = 0.02596
Variance sigma_squared = 3244.54
__________________________________________________________________
P(-508) = 0.02504
Variance sigma_squared = 3272.76
__________________________________________________________________
P(-512) = 0.0236
Variance sigma_squared = 3237.07
__________________________________________________________________
P(-516) = 0.02352
Variance sigma_squared = 3253.02
__________________________________________________________________
P(-520) = 0.02348
Variance sigma_squared = 3271.75
__________________________________________________________________
P(-524) = 0.02211
Variance sigma_squared = 3193.93
__________________________________________________________________
P(-528) = 0.022
Variance sigma_squared = 3185.49
__________________________________________________________________
P(-532) = 0.02076
Variance sigma_squared = 3257.23
__________________________________________________________________
P(-536) = 0.01927
Variance sigma_squared = 3197.19
__________________________________________________________________
P(-540) = 0.01808
Variance sigma_squared = 3268.26
__________________________________________________________________
P(-544) = 0.01742
Variance sigma_squared = 3242.76
__________________________________________________________________
P(-548) = 0.01629
Variance sigma_squared = 3245.94
__________________________________________________________________
P(-552) = 0.01467
Variance sigma_squared = 3188.52
__________________________________________________________________
P(-556) = 0.01429
Variance sigma_squared = 3239.2
__________________________________________________________________
P(-560) = 0.01301
Variance sigma_squared = 3229.09
__________________________________________________________________
P(-564) = 0.01249
Variance sigma_squared = 3251.9
__________________________________________________________________
P(-568) = 0.01226
Variance sigma_squared = 3202.81
__________________________________________________________________
P(-572) = 0.01045
Variance sigma_squared = 3235.04
__________________________________________________________________
P(-576) = 0.00998
Variance sigma_squared = 3255.13
__________________________________________________________________
P(-580) = 0.00902
Variance sigma_squared = 3271.82
__________________________________________________________________
P(-584) = 0.0081
Variance sigma_squared = 3246.92
__________________________________________________________________
P(-588) = 0.00821
Variance sigma_squared = 3270.59
__________________________________________________________________
P(-592) = 0.00779
Variance sigma_squared = 3265.64
__________________________________________________________________
P(-596) = 0.00616
Variance sigma_squared = 3229.64
__________________________________________________________________
P(-600) = 0.00589
Variance sigma_squared = 3172.27
__________________________________________________________________
P(-604) = 0.00568
Variance sigma_squared = 3242.89
__________________________________________________________________
P(-608) = 0.00472
Variance sigma_squared = 3250.42
__________________________________________________________________
P(-612) = 0.00366
Variance sigma_squared = 3272.59
__________________________________________________________________
P(-616) = 0.00377
Variance sigma_squared = 3277.35
__________________________________________________________________
P(-620) = 0.00329
Variance sigma_squared = 3244.88
__________________________________________________________________
P(-624) = 0.0025
Variance sigma_squared = 3237.85
__________________________________________________________________
P(-628) = 0.00241
Variance sigma_squared = 3292.5
__________________________________________________________________
P(-632) = 0.00194
Variance sigma_squared = 3199.75
__________________________________________________________________
P(-636) = 0.0014
Variance sigma_squared = 3142.61
__________________________________________________________________
P(-640) = 0.0015
Variance sigma_squared = 3210.86
__________________________________________________________________
P(-644) = 0.00122
Variance sigma_squared = 3238.5
__________________________________________________________________
P(-648) = 0.00078
Variance sigma_squared = 3148.36
__________________________________________________________________
P(-652) = 0.00105
Variance sigma_squared = 3293.37
__________________________________________________________________
P(-656) = 0.00069
Variance sigma_squared = 3242.9
__________________________________________________________________
P(-660) = 0.00056
Variance sigma_squared = 3254.05
__________________________________________________________________
P(-664) = 0.00051
Variance sigma_squared = 3244.99
__________________________________________________________________
P(-668) = 0.00045
Variance sigma_squared = 3232.27
__________________________________________________________________
P(-672) = 0.00038
Variance sigma_squared = 3276.1
__________________________________________________________________
P(-676) = 0.00015
Variance sigma_squared = 3176.92
__________________________________________________________________
P(-680) = 0.00013
Variance sigma_squared = 3151.96
__________________________________________________________________
P(-684) = 0.00015
Variance sigma_squared = 3243.91
__________________________________________________________________
P(-688) = 0.00022
Variance sigma_squared = 3334
__________________________________________________________________
P(-692) = 5e-05
Variance sigma_squared = 3147.06
__________________________________________________________________
P(-696) = 3e-05
Variance sigma_squared = 3285.78
__________________________________________________________________
P(-700) = 4e-05
Variance sigma_squared = 3229.82
__________________________________________________________________
247446 -493.769 0.99835
Computed variance = 3638.22
*/

    int idum;
    int n_spins, mcs, my_rank, numprocs;
    double average[5], total_average[5], initial_temp, final_temp, temp_step;


    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    //Lattice size, cycles
    n_spins = 100; mcs = 1000000; initial_temp = 2.20; final_temp = 2.40; temp_step = 0.005;

    MPI_Bcast (&n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int no_intervalls = mcs/numprocs;
    int myloop_begin = my_rank*no_intervalls + 1;
    int myloop_end = (my_rank+1)*no_intervalls;
    if ( (my_rank == numprocs-1) &&( myloop_end < mcs) ) myloop_end = mcs;

    //own seed to individual processors
    idum = -1-my_rank;
    srand(idum);
    // random starting point

    ofstream myfile;
    myfile.open("L100.txt");


    double T;
    double k = 1.38064852e-23;
    double J = 1;
    double E;
    double M;

    //loop over temperatures
    for (double That = initial_temp; That <= final_temp; That += temp_step)
    {
        //Initialize energy, magnetization and averages
        E = M = 0;
        T = abs(That*J/k);
        average[0] = 0; average[1] = 0; average[2] = 0; average[3] = 0; average[4] = 0;
        Ising L1 = Ising(J,n_spins,T);

        //monte carlo loop
        for (int cycles = myloop_begin; cycles <= myloop_end; cycles++)
        {
            L1.step_metropolis();
            E = L1.get_energy();
            M = L1.get_magnetization();
            average[0] += E; average[1] += E*E; average[2] += M; average[3] += M*M; average[4] += abs(M);
        }

        //find total average
        for( int i =0; i < 5; i++)
        {

            MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }

        if ( my_rank == 0)
        {

            myfile << That << " " << n_spins << " " << mcs << " " << total_average[0]/mcs << " " << total_average[1]/mcs << " "
                   << total_average[2]/mcs << " " << total_average[3]/mcs << " " << total_average[4]/mcs << endl;
        }

    }

    myfile.close();

    MPI_Finalize ();


    return 0;

}
