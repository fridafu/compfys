#include "planet.h"

//Here we make a planet
planet::planet(double mass, double x, double y, double z, double vx, double vy, double vz)
{
    this->mass = mass;
    position[0] = x;
    position[1] = y;
    position[2] = z;
    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;
    potential = 0.;
    kinetic = 0.;
}

//Find distance between two planets in 2D
double planet::distance(planet otherPlanet){ //Also making another planet
    double x1, x2, y1, y2, x_, y_;
    x1 = x; y1 = y;
    x2 = otherPlanet.position[0];
    y2 = otherPlanet.position[1];
    x_ = x2 - x1; y_ = y2 - y1;
    return sqrt(x_**2 + y_**2);
}

//Find the Gravitational force between two planets in 2D
double planet::GForce(planet otherPlanet, double G_constant){
    double r;
    r = distance(otherPlanet);
    return - (G_constant*(otherPlanet.mass)*mass)/(r**2);
}

//Find the acceleration of planet
double planet::Acceleration(planet otherPlanet, double G_constant){
    return GForce(otherPlanet, G_constant)/mass;
}

//Calculate Kinetic Energy of planet
double planet::KE(){
    v_squared = (velocity[0])**2 + (velocity[1])**2 + (velocity[2])**2;
    return 0.5*mass*v_squared;
}

//Caluclate Potential Energy of planet wrt other planet
double planet::PE(planet &otherPlanet, double G_const, double epsilon){
    return mass*G_const*distance(otherPlanet);
}
