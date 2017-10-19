#ifndef SOLVER_H
#define SOLVER_H
#include "planet.h"
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include "armadillo"
#define _USE_MATH_DEFINES
//#define pi M_PI
using std::vector;
using namespace arma;

class solver
{
public:
    friend class planet;

    // what kind of properties are there in the solver
    //double radius, mass_tot, G;
    //int Nplanets; // number of planets in the system we look at

    // the energy in the system
    int planets_tot = 0;
    double radius = 100;
    double mass_tot = 0;
    double G = 4*M_PI*M_PI;
    double kinetic_tot = 0;
    double potential_tot = 0;

    // initializers
    solver();
    solver(double radi);

    // functions in the class solver
    //void add(planet newplanet);
    void G_constant();
    //void print_trajectory();//std::ofstream &output, int dimensions, double time, int number);
    void velVerlet(int dimensions, double h, int Nplanets);



};

#endif // SOLVER_H
