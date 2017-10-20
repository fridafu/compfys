#ifndef SOLVER_H
#define SOLVER_H

#include <planet.h>
#include <vector>
#include <armadillo>



using namespace arma;
using namespace std;

class Solver
{
private:
    double pi;
    std::vector<planet> objects;
    int numobj;
    vec t;
    double dt;
    vec accel2;
    bool sunfixed;
    mat accel1;
    ofstream myfile;
    ofstream mytimes;

public:
    Solver();
    void addPlanet(planet newplanet);
    void solve(vec times, bool sunfix); //set sunfix = true if last planet added should be held fixed.
    void stepVerlet();

};

#endif // SOLVER_H
