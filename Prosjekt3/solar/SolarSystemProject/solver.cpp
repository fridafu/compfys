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


void Solver::solve(vec times, bool sunfix, bool writefile, int skipwrite, int method)
// times - vector of times.
//set sunfix = true if the last planet added should be fixed in position
// set writefile = true if you want the positions written to file
// method = 1 for verlet. method = 2 for euler
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
            myfile << objects[i].position(0) << " " << objects[i].position(1) << " " << objects[i].position(2) << " ";
        }
        myfile << endl;
    }
    for (long int i = 0; i < size(t)(0); i++)
    {
        if (method == 1)
        {
            stepVerlet();
        }
        if (method == 2)
        {
            stepEuler();
        }

        if (writefile)
        {
            if (i%skipwrite == 0)
            {
                mytimes << t(i) << endl;
                cout << 100*float(i)/size(t)(0) << endl;

                for (int j = 0; j < numobj; j++)
                {
                    myfile << objects[j].position(0) << " " << objects[j].position(1) << " " << objects[j].position(2) << " ";

                }
                myfile << endl;
            }

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

void Solver::testVel(vec times, bool sunfix, vec initvel, int method)
{

    sunfixed = sunfix;
    if (sunfixed)
    {
        numobj--;
    }

    double tol = 0.03;
    long int i;
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

        t = times;
        dt = t(1) - t(0);

            for (long int j = 0; j < size(t)(0); j++)
            {
                if (method == 1)
                {
                    stepVerlet();
                }
                if (method == 2)
                {
                    stepEuler();
                }
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


void Solver::testStability(bool sunfix, vec dt_, int method)
{

    sunfixed = sunfix;
    if (sunfixed)
    {
        numobj--;
    }


    mydt.open("dt.txt");
    mydistance.open("distance.txt");

    double x = objects[0].position(0);
    double y = objects[0].position(1);
    double z = objects[0].position(2);
    double vx = objects[0].velocity(0);
    double vy = objects[0].velocity(1);
    double vz = objects[0].velocity(2);


    for (int i = 0; i < size(dt_)(0); i++)
    {


        set_dt(dt_(i));
        mydt << dt_(i) << endl;
        objects[0].position(0) = x;
        objects[0].position(1) = y;
        objects[0].position(2) = z;
        objects[0].velocity(0) = vx;
        objects[0].velocity(1) = vy;
        objects[0].velocity(2) = vz;




        long int steps = 1000/dt_(i);

        firststep = true;

        for (long int k = 0; k < steps; k++)
        {

            if (method == 1)
            {
                stepVerlet();
            }
            if (method == 2)
            {
                stepEuler();
            }


            if (k%1000000 == 0)
            {
                cout << 100*float(k)/steps << endl;
            }
        }
        mydistance << objects[0].distance(objects[1]) << endl;


    }
    mydt.close();
    mydistance.close();
}

void Solver::testConservation(bool sunfix, int method)
{
    sunfixed = sunfix;
    if (sunfixed)
    {
        numobj--;
    }

    vec t = linspace(0,10,10000);
    dt = t(1) - t(0);

    double KE_end;
    double PE_end;
    double L_end;
    double tol = 0.0001;
    vec KE_ = zeros(size(t));
    vec PE_ = zeros(size(t));
    vec E_ = zeros(size(t));
    vec L_ = zeros(size(t));

    int i;
    for(i=0;i<10000;i++){
        //objects[0].position(0) = 1.;
        //objects[0].position(1) = 0.0;
        //objects[0].position(2) = 0.0;
        //objects[0].velocity(0) = 0.0;
        //objects[0].velocity(2) = 0.0;
        //objects[0].velocity(1) = 6.3;

        if (method == 1)
        {
            stepVerlet();
        }
        if (method == 2)
        {
            stepEuler();
        }
        KE_(i) = objects[0].KE();
        PE_(i) = objects[0].PE(objects[1], 4*pi*pi);
        E_(i) = KE_(i) + PE_(i);
        L_(i) = objects[0].L(objects[1]);
        }
    KE_end = KE_.end()[-2];
    PE_end = PE_.end()[-2];
    L_end = L_.end()[-2];
    cout << KE_end - tol << KE_(1) << KE_end + tol << endl;
    if(KE_end - tol < KE_(1)  && KE_(1) < KE_end + tol){
           cout << "Kinetic energy is conserved! (Circular orbit!)" << endl;
    }

    if(PE_end - tol < PE_(1) && PE_(1) < PE_end + tol){
           cout << "Potential energy is conserved! (Circular orbit!)" << endl;
    }

    if((PE_end + KE_end) - tol < (PE_(1) + KE_(1)) && (PE_(1) + KE_(1)) < (PE_end + KE_end) + tol){
           cout << "Total energy is conserved!" << endl;
    }

    if(L_end - tol < L_(1) && L_(1) < L_end + tol){
           cout << "Angular momentum is conserved!" << endl;
    }
}
