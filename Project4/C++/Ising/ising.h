#ifndef ISING_H
#define ISING_H

#include "armadillo"

using namespace arma;
using namespace std;


class Ising
{
public:
    Ising(double coupling, int l, double temp);
    //coupling is a constant expressing the strength of the interaction between neighboring spins
    //l is the dimension of the 2D lattice (lxl)
    //temp is the temperature in kelvin.
    //When initializing object, a random state is generated.

    mat rand_state();
    //generates and returns a random state

    void set_state(mat S);
    //takes a matrix S as input and sets it as the state.

    void energy();
    //calculates the energy and stores it in variable E

    mat flip_rand_spin(mat S);
    //flips random spin in input matrix.

    void step_metropolis();
    //performs one step of the metropolis algorithm
    //stores the number of times this function is used in variable stepcount
    //stores the number of accepted configurations in variable totaccept

    void exp_vals(int steps);
    //calculates expectation values with 'steps' steps.

    void magnetization();
    //calculates the magnetization and stores it in variable M

    double heat_capacity();
    //returns the heat capacity

    double magnetic_susceptibility();
    //returns the magnetic susceptibility

    vec get_expectation_values();
    //returns expectation values {mean energy, mean energy squared, mean magnetization, mean magnetization squared, mean absolute magnetization}

    double get_energy();
    //returns the energy of the system

    int get_configurations();
    //returns the total accepted configurations.

    void reset_expectation_values();
    //resets the expectation values

    double energy_probability(int inital_steps, int steps, int En);
    //This function will do 'inital_steps' steps with the metropolis algorithm, then count the number of times energy En occurs
    //and return the probability P(En).

    double get_magnetization();
    //returns the magnetization of the system


private:

    int totaccept;
    double k;
    double r;
    double E;
    double M;
    long double expE;
    long double expE2;
    long double expM;
    long double expM2;
    long double expabsM;
    double dM;
    int up;
    int down;
    int right;
    int left;
    int fr;
    int fc;
    double T;
    double beta;
    mat state;
    int L;
    double J;
    double dE;
    long int stepcount;
};

#endif // ISING_H
