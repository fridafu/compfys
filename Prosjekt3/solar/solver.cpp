#include "solver.h"



solver::solver()
{
    int planets_tot = 0;
    double radius = 100;
    double mass_tot = 0;
    double G = 4*M_PI*M_PI;
    double kinetic_tot = 0;
    double potential_tot = 0;
}

solver::solver(double radi){
    planets_tot = 0;
    radius = radi;
    mass_tot = 0;
    G = 4*M_PI*M_PI;

}

double solver::G_constant(){ //double mass4pi){
    return 4*M_PI*M_PI;// mass4pi; //radius*radius*radius;

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

void solver::velVerlet(double h, int n, planet* p1, planet* p2){// call on this in a for loop in overlooking class
    vec a(3);
    vec anew(3);
    // Velocity Verlet method for trajectory calculations
    a = p1->Acceleration(*p2,G_constant())%p1->position;
    p1->position += h*p1->velocity + 0.5*h*h*a;
    anew = p1->Acceleration(*p2,G_constant())%p1->position;
    p1->velocity += 0.5*h*(a+anew);
}

void solver::add(planet newplanet)
{
    total_planets += 1;
    total_mass += newplanet.mass;
    all_planets.push_back(newplanet);
}

void solver::addM(planet newplanet)
{
    total_planets +=1;
    all_planets.push_back(newplanet);
}


/*
axi = -GM_sun*x(i)/pow(rsqrt,1.5);
ayi = -GM_sun*y(i)/pow(rsqrt,1.5);
x(i+1) = x(i) + h*vx(i) + 0.5*h*h*axi;
y(i+1) = y(i) + h*vy(i) + 0.5*h*h*ayi;
axi1 = -GM_sun*x(i+1)/pow(rsqrt,1.5); // calculate ax(i+1) and
ayi1 = -GM_sun*y(i+1)/pow(rsqrt,1.5); // calculate ax(i+1) and
vx(i+1) = vx(i) + 0.5*h*(axi+axi1);
vy(i+1) = vy(i) + 0.5*h*(ayi+ayi1);
//myfile << x(i+1) << "  " << y(i+1) << endl; // write rest of x, y coordinates to file
}
myfile.close();
*/

/*
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
*/
