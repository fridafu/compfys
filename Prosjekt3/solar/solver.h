#ifndef SOLVER_H
#define SOLVER_H
#include "planet.h"
#include <vector>
#include <fstream>
using std::vector;

class solver
{
public:
    friend class planet;

    // what kind of properties are there in the solver
    double radius, mass_tot, G;
    int Nplanets; // number of planets in the system we look at


    // the energy in the system
    double kinetic_tot, potential_tot;


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
