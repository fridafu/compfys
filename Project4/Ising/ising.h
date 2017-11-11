#ifndef ISING_H
#define ISING_H

#include "armadillo"

using namespace arma;
using namespace std;


class Ising
{
public:
    Ising(double coupling, int l, double temp);
    mat rand_state();
    void set_state(mat S);
    void energy();
    mat flip_rand_spin(mat S);
    void step_metropolis();
    void exp_vals(int steps);
    void magnetization();
    double heat_capacity();
    double magnetic_susceptibility();
    vec get_expectation_values();
    double get_energy();
    int get_configurations();
    void reset_expectation_values();
    double energy_probability(int inital_steps, int steps, int En);


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
