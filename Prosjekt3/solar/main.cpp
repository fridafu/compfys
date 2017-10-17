#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#define pi M_PI
#include "solver.h"
#include "planet.h"

using namespace std;
using namespace arma;

int forward_Euler(int n, double h, vec x, vec y, vec vx, vec vy, string filename, double GM_sun=4*pi*pi);
int velocity_Verlet(int n, double h, vec x, vec y, vec vx, vec vy, string filename, double GM_sun=4*pi*pi);

int main()
{
    double x0 = 1.0; // [AU]
    double y0 = 0.0; // [AU]
    double t_max = 3.0; // [yrs]
    double t_min = 0.0; // [yrs]
    int N = 300000; // number of steps
    double h = (t_max - t_min)/(N); // Step size h
    // initial conditions for postions
    vec x = zeros<vec>(N);
    vec y = zeros<vec>(N);
    vec vx = zeros<vec>(N);
    vec vy = zeros<vec>(N);
    x(0) = x0;
    y(0) = y0;
    string filename1, filename2;
    cout << "Forward Euler. Give me a filname: fe... .txt    ";
    cin >> filename1;
    cout << "Velocity Verlet. Give me a filname: vv... .txt    ";
    cin >> filename2;

    vx(0) = 0.0;//-GM_sun*x0/pow(r0,1.5);
    vy(0) = 2*pi;//GM_sun*y0/pow(r0,1.5);

    // calculate the trajectories using both methods
    forward_Euler(N, h, x, y, vx, vy, filename1);
    velocity_Verlet(N, h, x, y, vx, vy, filename2);
    return 0;
}

int forward_Euler(int n, double h, vec x, vec y, vec vx, vec vy, string filename, double GM_sun){
    ofstream myfile;
    myfile.open("/home/hannahcb/compfys/Prosjekt3/Python/fe" +  filename +".txt");
    myfile << x(0) << "  " << y(0) << endl; // write initial conditions to file
    double ax, ay, rsqrt; // rsqrt = (x^2 + y^2)^0.5
    // Forward Euler to calculate trajectory
    for (int i = 0; i<n-1;i++){
        rsqrt = (x(i)*x(i)+y(i)*y(i));
        ax = -GM_sun*x(i)/pow(rsqrt,1.5);
        ay = -GM_sun*y(i)/pow(rsqrt,1.5);
        x(i+1) = x(i) + h*vx(i);
        y(i+1) = y(i) + h*vy(i);
        vx(i+1) = vx(i) + h*ax;
        vy(i+1) = vy(i) + h*ay;
        myfile << x(i+1) << "  " << y(i+1) << endl; // write rest of x, y coordinates to file
    }
    myfile.close();
    return 0;
}

int velocity_Verlet(int n, double h, vec x, vec y, vec vx, vec vy, string filename, double GM_sun){
    ofstream myfile;
    myfile.open("/home/hannahcb/compfys/Prosjekt3/Python/vv" +  filename +".txt");
    myfile << x(0) << "  " << y(0) << endl; // write initial conditions to file
    double axi, axi1, ayi, ayi1, rsqrt;
    for (int i = 0; i<n-1;i++){ // Velocity Verlet method for trajectory calculations
        rsqrt = (x(i)*x(i)+y(i)*y(i));
        axi = -GM_sun*x(i)/pow(rsqrt,1.5);
        ayi = -GM_sun*y(i)/pow(rsqrt,1.5);
        x(i+1) = x(i) + h*vx(i) + 0.5*h*h*axi;
        y(i+1) = y(i) + h*vy(i) + 0.5*h*h*ayi;
        axi1 = -GM_sun*x(i+1)/pow(rsqrt,1.5); // calculate ax(i+1) and
        ayi1 = -GM_sun*y(i+1)/pow(rsqrt,1.5); // calculate ax(i+1) and
        vx(i+1) = vx(i) + 0.5*h*(axi+axi1);
        vy(i+1) = vy(i) + 0.5*h*(ayi+ayi1);
        myfile << x(i+1) << "  " << y(i+1) << endl; // write rest of x, y coordinates to file
    }
    myfile.close();
    return 0;
}

