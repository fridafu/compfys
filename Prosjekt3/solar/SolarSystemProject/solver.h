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
    mat accel2;
    bool sunfixed;
    mat accel1;
    ofstream myfile;
    ofstream mytimes;
    ofstream mydt;
    ofstream mydistance;
    double circvel;
    bool firststep;

public:
    Solver();
    void addPlanet(planet newplanet);
    void solve(vec times, bool sunfix, bool writefile = true, int skipwrite = 1, string method = "1"); //set sunfix = true if last planet added should be held fixed.
    void stepVerlet();
    void stepEuler();
    void testVel(vec times, bool sunfix, vec initvel);
    void testStability(bool sunfix, vec dt_);
    void set_dt(double deltat);
    void set_sunfixed(bool sunfix);
    vec get_position(int planetnr);
    void testConservation(bool sunfix);

};

#endif // SOLVER_H
