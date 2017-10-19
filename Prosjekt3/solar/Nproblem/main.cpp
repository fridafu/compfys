#include <iostream>
#include "planet.h"
#include <fstream>
#include "armadillo"
using namespace std;
using namespace arma;


template <class T>

class Nproblem {
    private:
        T objects[];
        int numobj;
        vec t;
        double dt;
        vec accel1(numobj);
        double accel2;
    public:

        Nproblem(T obj[]) {
            objects = obj;
            numobj = sizeof(obj);




        }
        double solve(vec times) {
            t = times;
            dt = t(1) - t(0);
            ofstream mytimes;
            mytimes.open("times.txt");
            for (int i = 0; i < sizeof(t); i++) {
                mytimes << t(i) << endl;
            }
            mytimes.close();
            this->advance();
        }

}

class Verlet : public Nproblem {
    public:
        void advance() {

           ofstream myfile;
           myfile.open ("solarsystem.txt")

           for (int i = 0; i < sizeof(t); i++) {
               for (int j = 0; j < numobj - 1; j++) {

                   vec force = zeros(3);

                   for (int k = 0; k < numobj - 1; k++) {
                       force += objects(j).Force(objects(k));
                   }
                   force += objects(j).Force(objects(numobj-1));

                   accel1(j) = force/objects(j).mass;
                   objects(j).position = objects(j).position + objects(j).velocity*dt + 0.5*accel1(j)*dt**2;
                   myfile << objects(j).position + " ";
               }
               for (int j = 0; j < numobj - 1; j++) {
                   force = zeros(3);


                   for (int k = 0; k < numobj - 1; k++) {
                       force += objects(j).Force(objects(k));
                   }
                   force += objects(j).Force(objects(numobj - 1));

                   accel2 = force/objects(j).mass;
                   objects(j).velocity = objects(j).velocity + 0.5*(accel1(j) + accel2)*dt;
               }
               myfile << endl;
            }
            myfile.close();

}

int main()
{
    solarsystem = listeavplanetobjekter;
    solver = Verlet(solarsystem);
    solver.solve();

    return 0;
}
