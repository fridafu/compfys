#include <iostream>
#include <armadillo>
#include <planet.h>
#include <solver.h>


using namespace std;
using namespace arma;



int main()
{

    //3g - Relativistic correction - Angle

    planet mercury(0.00044, -0.3075, 0, 0, 0, -12.44, 0);
    planet sun(1, 0, 0, 0, 0, 0, 0);

    mercury.relcheck = true; //Includes relativistic term in force.

    Solver mercurysun;

    mercurysun.addPlanet(mercury);

    mercurysun.addPlanet(sun);

    ofstream relfile;
    ofstream reltime;

    relfile.open("relfile2.txt");
    reltime.open("reltime2.txt");

    double dt = 1E-8;
    double T = 100;
    mercurysun.set_dt(dt);
    mercurysun.set_sunfixed(true);
    int u = 1;
    double r1;
    double r2;
    double r3;
    vec pos2;

    for (long int i = 0; i < T/(dt); i++)
    {
        mercurysun.stepVerlet();
        r1 = norm(mercurysun.get_position(0));

        if ( (r2 < r3) && (r2 < r1) && (i > 3) )
        {
            cout << "orbit " << u << endl;
            relfile << atan(pos2(1)/pos2(0))*206264.806 << endl;
            reltime << i*dt << endl;
            u++;
        }

        r3 = r2;
        r2 = r1;
        pos2 = mercurysun.get_position(0);
    }
    relfile.close();
    reltime.close();




}
