#include "planet.h"

//Here we make a planet
planet::planet(double mass, double x, double y, double z, double vx, double vy, double vz)
{
    this->mass = mass;
    position(0) = x;
    position(1) = y;
    position(2) = z;
    velocity(0) = vx;
    velocity(1) = vy;
    velocity(2) = vz;
    potential = 0.;
    kinetic = 0.;
}

//Find distance between two planets in 2D
double planet::distance(planet otherPlanet){ //Also making another planet
    vec pos1,pos2;
    double r;

    pos1 = this->position;
    pos2 = otherPlanet.position;
    r = norm(pos1-pos2);
    return r;
 }



//Find the Gravitational force between two planets in 2D
vec planet::GForce(planet otherPlanet, double G_constant){
    double r;
    r = distance(otherPlanet);
    return -(position - otherPlanet.position)*(G_constant*(otherPlanet.mass)*this->mass)/(r*r*r);
}

//Find the acceleration of planet
vec planet::Acceleration(planet otherPlanet, double G_constant){
    return GForce(otherPlanet, G_constant)/this->mass;
}

//Calculate Kinetic Energy of planet
double planet::KE(){
    return 0.5*mass*dot(velocity,velocity);
}

//Caluclate Potential Energy of planet wrt other planet
double planet::PE(planet &otherPlanet, double G_const, double epsilon){
    return mass*G_const*distance(otherPlanet);
}


