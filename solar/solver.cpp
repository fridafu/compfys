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
    kinetic_tot = 0;
    potential_tot = 0;
}

void solver::G_constant(){
    G = (4*pi*pi/32)*radius*radius*radius/mass_tot;
}

void solver::print_trajectory(std:ofstram &output, string filename, int N, int dimensions, double time){ // write mass, trajectory and velocity to file 'output'
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


void velVerlet(){


    // define time step size h = t_final - t_initial / integration points
    double h = final_time/((double) N);
    double time = 0.0;


    // set up arrays for positions





}

}
