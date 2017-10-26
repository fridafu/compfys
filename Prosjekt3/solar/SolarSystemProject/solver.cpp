#include "solver.h"

Solver::Solver()
{
    pi = datum::pi;
    numobj = 0;
    firststep = true;

}
void Solver::addPlanet(planet newplanet)
//adds a planet object.
{
    numobj++;
    objects.push_back(newplanet);
    mat accel1_temp(3,numobj);
    accel1 = accel1_temp;
    mat accel2_temp(3,numobj);
    accel2 = accel2_temp;
}

void Solver::set_dt(double deltat)
//Set dt when manually using stepVelvet or stepEuler - function
{
    dt = deltat;
}
void Solver::set_sunfixed(bool sunfix)
//Use this to set last planet fixed when manually using stepVelvet or stepEuler - function.
{
    sunfixed = sunfix;
    if (sunfixed)
    {
        numobj--;
    }

}

vec Solver::get_position(int planetnr)
//return position of planet number "planetnr"
{
    return objects[planetnr].position;
}

void Solver::solve(vec times, bool sunfix, bool writefile, string method)

// times - vector of times.
//set sunfix = true if the last planet added should be fixed in position
// set writefile = true if you want the positions written to file
// method = "1" for verlet. method = "2" for euler

{

    sunfixed = sunfix;
    if (sunfixed)
    {
        numobj--;
    }


    t = times;
    dt = t(1) - t(0);

    mytimes.open("times.txt");
    myfile.open ("solarsystem.txt");

    if (writefile)
    {
        for (int i = 0; i < numobj; i++)
        {
            myfile << objects[i].position(0) << " " << objects[i].position(1) << " " << objects[i].position(2) << " ";;
        }
        myfile << endl;
    }
    for (int i = 0; i < size(t)(0); i++)
    {
        if (method == "1")
        {
            stepVerlet();
        }
        if (method == "2")
        {
            stepEuler();
        }

        if (writefile)
        {
            mytimes << t(i) << endl;


            for (int j = 0; j < numobj; j++)
            {
                myfile << objects[j].position(0) << " " << objects[j].position(1) << " " << objects[j].position(2) << " ";
            }
            myfile << endl;
        }

    }


    mytimes.close();
    myfile.close();

}



void Solver::stepVerlet()
//Step forward using velocity verlet
{
    if (firststep)
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

        }
        firststep = false;
    }
    for (int j = 0; j < numobj; j++)
    {
        objects[j].position = objects[j].position + objects[j].velocity*dt + 0.5*accel1.col(j)*dt*dt;
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

        accel2.col(j) = force/objects[j].mass;


    }
    for (int j = 0; j < numobj; j++)
    {
        objects[j].velocity = objects[j].velocity + 0.5*(accel1.col(j) + accel2.col(j))*dt;
    }
    accel1 = accel2;
}

void Solver::stepEuler()
//Step forward using euler method.
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


    }
    for (int j = 0; j < numobj; j++)
    {
        objects[j].position = objects[j].position + objects[j].velocity*dt;
        objects[j].velocity = objects[j].velocity + accel1.col(j)*dt;
    }
}

void Solver::testVel(vec times, bool sunfix, vec initvel)
{
    double tol = 0.03;
    int i;
    for(i=0; i<30; i++){
        int bigv = 0;
        double vel;
        vel = initvel(i);
        objects[0].position(0) = 1.;
        objects[0].position(1) = 0.0;
        objects[0].position(2) = 0.0;
        objects[0].velocity(0) = 0.0;
        objects[0].velocity(2) = 0.0;
        objects[0].velocity(1) = vel;

        sunfixed = sunfix;
        if (sunfixed)
        {
            numobj--;
        }

        t = times;
        dt = t(1) - t(0);

            for (int j = 0; j < size(t)(0); j++)
            {
                stepVerlet();
                double r_ = abs(objects[0].distance(objects[1]));
                if(1.- tol > r_ or r_ > 1.+tol){
                       bigv += 1;
                }
            }
            if(bigv == 0){
                circvel = vel;
            }
        }
    cout << "The initial speed for a circular orbit is = " << circvel << " AU/years" << endl;
}


void Solver::testStability(bool sunfix){
    sunfixed = sunfix;
    if (sunfixed)
    {
        numobj--;
    }

    int i;
    double n;
    vec tpoints = linspace(100,1000000,100);

    vec distances = zeros(500);
    vec dt_ = linspace(1e-7,0.5,500);
    mydt.open("dt.txt");
    mydistance.open("distance.txt");
    for(i=0;i<500;i++){
        objects[0].position(0) = 1.;
        objects[0].position(1) = 0.0;
        objects[0].position(2) = 0.0;
        objects[0].velocity(0) = 0.0;
        objects[0].velocity(2) = 0.0;
        objects[0].velocity(1) = 6.3;

        //n = tpoints(i);
        dt = dt_(i);
        mydt << dt << endl;;
        vec r_ = zeros(500);
        for (int k = 0; k < size(dt_)(0); k++){
            stepVerlet();
            r_(k) = objects[0].distance(objects[1]);
            }
        distances(i) = r_.end()[-2];
        mydistance << distances(i) << endl;;
        }
    mydt.close();
    mydistance.close();
}


