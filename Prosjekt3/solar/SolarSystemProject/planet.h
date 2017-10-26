#ifndef PLANET_H
#define PLANET_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
using std::vector;
#include "armadillo"
using namespace arma;

class planet
{
public:
    // Properties
    double mass;
    vec position = zeros(3);
    vec velocity = zeros(3);
    vec pos1,pos2;
    double potential;
    double kinetic;
    double beta;
    bool relcheck;

    // Initializers
    planet();
    planet(double mass,double x,double y,double z,double vx, double vy,double vz);

    // Functions
    double distance(planet& otherPlanet);
    vec GForce(planet& otherPlanet, double G_constant);
    vec Acceleration(planet& otherPlanet, double G_constant);
    double KE();



    double PE(planet& otherPlanet, double G_constant);

    double L(planet& otherPlanet);



};

#endif // PLANET_H
