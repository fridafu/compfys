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
    double energy();
    mat flip_rand_spin(mat S);
    void step_exp_vals();
    vec exp_vals(int steps);
    double magnetization();
    double heat_capacity();
    double magnetic_susceptibility();
    double E;
    double M;
    double E2;
    double M2;
    double absM;

private:
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
    double sumdE;
    double sumdM;
    mat state;
    int L;
    double J;

    double dE;
    int stepcount;
};

#endif // ISING_H
