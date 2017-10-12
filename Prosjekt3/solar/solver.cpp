#include "solver.h"
#include <armadillo>
#include <cmath>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#define pi M_PI

solver::solver()
{
    planets_tot = 0;
    radius = 100;
    mass_tot = 0;
    G = 4*pi*pi;
    kinetic_tot = 0;
    potential_tot = 0;
}

solver::solver(double radi){
    planets_tot = 0;
    radius = radi;
    mass_tot = 0;
    G = 4*pi*pi;

}

void solver::G_constant(){ //double mass4pi){
    G = 4*pi*pi// mass4pi; //radius*radius*radius;
}

/*void solver::print_trajectory(std:ofstram &output, string filename, int N, int dimensions, double time){ // write mass, trajectory and velocity to file 'output'
    if (dimensions > 3 || dimensions <= 0) dimensions = 3;

    else{
        for (int i = 0; i < N-1; i++){
            planet &current = all_planets[i];
            output << time << "\t" << i+1 << "\t" << current.mass;
            for (int = j; j < dimensions; j++) output << "\t" << current.positions[j];
            for (int = j; j < dimensions; j++) output << "\t" << current.velocity[j];
            output << std::endl;
        }
    }

*/
void solver::velVerlet(int dimensions, double h, int Nplanets){


    // define time step size h = t_final - t_initial / integration points
    //double h = final_time/((double) N);
    //double time = 0.0;

    // make arrays for position,

    vec position = zeros<vec>(dimensions);
    vec velocity = zeros<vec>(dimensions);
    vec acceleration = zeros<vec>(dimensions);// need to initialize the start
    vec newacceleration = zeros<vec>(dimensions);
    // planets = [Sun, Earth,Mars]
    for(int j=0;j<Nplanets;j++){
        planets(j).position += h*planets(j).velocity + 0.5*h*h*planets(j).acceleration;
        planets(j).newacceleration = zeros<vec>(dimensions); //somefunctionthatcalculatesacceleration();
        planets(j).velocity += 0.5*h*(planets(j).acceleration+planets(j).newacceleration);
        planets(j).acceleraton = planets(j).newacceleration;

    }
    for(int i = 0; i<N-1; i++){




        // call on print function
        // print to file?

    }



    // set up arrays for positions

}
