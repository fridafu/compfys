#include <iostream>
#include "planet.h"
using namespace std;


template <class T>

class Nproblem {
    private:
        T objects[];
        int numobj;
        double t[];
        double rprev[numobj][3];
        double masses[numobj];
        double v0[numobj][3];
    public:
        Nproblem(T obj[]) {
            numobj = sizeof(obj);
            objects = obj;

            for (int n = 0; n < numobj; n++) {
                rprev[n] = objects(n).position;
                masses[n] = objects(n).mass;


            }
        }
        double solve(double times[]) {
            double r = new double[sizeof(times)][numobj][3];
            r[0] = rprev;

            t = times[];
            for (int i = 1; i < sizeof(t); i++) {
                for (int j = 0; j < numobj; j++) {
                    double Fnext[3] = {0,0,0};
                    for (int k = 0; k < numobj; k++) {
                        if (k != j) {
                            Fnext += objects(j).Gforce(objects(k));

                        }
                    acceleration = Fnext/mass(j);
                    r[i][j] = advance();

                    }
                for (int l = 0; l < numobj; l++) {
                    object(l).set_position(r[i][l])
                }


                }
                rprev = r[i]
            }


        }

}

class Verlet : public Nproblem {
    public:
        double[] advance(double[] acceleration) {


        }
}

int main()
{
    cout << "Hello World!" << endl;
    return 0;
}
