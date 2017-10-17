#ifndef PLANET_H
#define PLANET_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
using std::vector;


class planet
{
public:
    // Properties
    double mass;
    double position[3];
    double velocity[3];
    double potential;
    double kinetic;

    // Initializers
    planet();
    planet(double mass,double x,double y,double z,double vx, double vy,double vz);

    // Functions
    double distance(planet otherPlanet);
    double GForce(planet otherPlanet, double G_constant);
    double Acceleration(planet otherPlanet, double G_constant);
    double KE();
    double PE(planet &otherPlanet, double G_constant, double epsilon);

};

#endif // PLANET_H
