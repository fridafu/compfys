#include <iostream>
#include "planet.h"
#include "Nproblem.h"
#include <fstream>
#include "armadillo"
using namespace std;
using namespace arma;




Nproblem::Nproblem(planet obj)
{
    objects = obj;
    numobj = sizeof(obj);
    T objects[];
    int numobj;
    vec t;
    double dt;
    vec accel1(numobj);
    double accel2;
}
double Nproblem::solve(vec times)
{
    t = times;
    dt = t(1) - t(0);
    ofstream mytimes;
    mytimes.open("times.txt");
    ofstream myfile;
    myfile.open ("solarsystem.txt");
    myfile << objects(j).position + " ";
    for (int i = 0; i < sizeof(t); i++)
    {
        mytimes << t(i) << endl;
        this->stepVerlet();

    }
    mytimes.close();
    myfile.close();

}


{
void Nproblem::stepVerlet()
{


        for (int j = 0; j < numobj - 1; j++)
        {

        vec force = zeros(3);

        for (int k = 0; k < numobj - 1; k++)
        {
            if (j != k)
            {
            force += objects(j).GForce(objects(k), 4*pi*pi);
            }
        }
        force += objects(j).GForce(objects(numobj-1), 4*pi*pi);

        accel1(j) = force/objects(j).mass;
        objects(j).position = objects(j).position + objects(j).velocity*dt + 0.5*accel1(j)*dt*dt;
        myfile << objects(j).position + " ";
        }
        for (int j = 0; j < numobj - 1; j++)
        {
            force = zeros(3);


            for (int k = 0; k < numobj - 1; k++)
            {
                if (j != k)
                {
                force += objects(j).GForce(objects(k), 4*pi*pi);
                }
            }
            force += objects(j).GForce(objects(numobj - 1), 4*pi*pi);

            accel2 = force/objects(j).mass;
            objects(j).velocity = objects(j).velocity + 0.5*(accel1(j) + accel2)*dt;
        }
        myfile << endl;
}


