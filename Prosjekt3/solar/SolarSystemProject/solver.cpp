#include "solver.h"




Solver::Solver()
{
    pi = datum::pi;
    numobj = 0;
    vec accel2_temp(3);
    accel2 = accel2_temp;
}
void Solver::addPlanet(planet newplanet)
{
    numobj++;
    objects.push_back(newplanet);
    mat accel1_temp(3,numobj);
    accel1 = accel1_temp;
}

void Solver::solve(vec times, bool sunfix)
{
    sunfixed = sunfix;
    if (sunfixed)
    {
        numobj--;
    }

    string method = "1";
    t = times;
    dt = t(1) - t(0);

    mytimes.open("times.txt");
    myfile.open ("solarsystem.txt");

    for (int i = 0; i < numobj; i++)
    {
        myfile << objects[i].position(0) << " " << objects[i].position(1) << " " << objects[i].position(2) << " ";;
    }
    myfile << endl;
    if (method == "1")
    {
        for (int i = 0; i < size(t)(0); i++)
        {
            mytimes << t(i) << endl;

            stepVerlet();

        }
    }
    mytimes.close();
    myfile.close();

}



void Solver::stepVerlet()
{

    for (int j = 0; j < numobj; j++)
    {

        vec force = zeros(3);

        for (int k = 0; k < numobj; k++)
        {
            if (j != k)
            {
            force += objects[j].GForce(objects[k], 4*pi*pi);
            }
        }
        if (sunfixed)
        {
            force += objects[j].GForce(objects[numobj], 4*pi*pi);
        }

        accel1.col(j) = force/objects[j].mass;
        objects[j].position = objects[j].position + objects[j].velocity*dt + 0.5*accel1.col(j)*dt*dt;

        myfile << objects[j].position(0) << " " << objects[j].position(1) << " " << objects[j].position(2) << " ";

    }
    for (int j = 0; j < numobj; j++)
    {
        vec force = zeros(3);

        for (int k = 0; k < numobj; k++)
        {
            if (j != k)
            {
            force += objects[j].GForce(objects[k], 4*pi*pi);
            }
        }
        if (sunfixed)
        {
           force += objects[j].GForce(objects[numobj], 4*pi*pi);
        }

        accel2 = force/objects[j].mass;

        objects[j].velocity = objects[j].velocity + 0.5*(accel1.col(j) + accel2)*dt;
    }
    myfile << endl;
}


