#include <iostream>
#include <armadillo>
#include <planet.h>
#include <solver.h>

using namespace std;
using namespace arma;

int main()
{
    double diy = 365.242199;


    planet earth(0.000003,8.930309961463524E-01, 4.508411860073833E-01, -1.549928799609221E-04, -diy*7.978069853256020E-03, diy*1.533806773162681E-02, -diy*7.469966577773285E-07); // Earth: (mass,x,y,z,vx,vy,vz)
    planet mars(3.3E-7, -1.575528590018643E+00, 5.364921619154235E-01, 4.971260736828571E-02, -diy*3.942819907858700E-03, -diy*1.206234082100237E-02, -diy*1.561289187128983E-04);
    planet saturn(0.000275, -3.369277467674985E-01, -1.004985402690239E+01, 1.881475188874106E-01, diy*5.269910272418428E-03, -diy*2.053763749248873E-04, -diy*2.059570456178141E-04);
    planet uranus(0.00044, 1.785433444765009E+01, 8.819856193817925E+00, -1.985484452135033E-01, -diy*1.770719891211758E-03, diy*3.342961457662236E-03, diy*3.539795634394496E-05);
    planet jupiter(0.00095, -4.572777635194016E+00, -2.939093897020645E+00, 1.144702114470432E-01, diy*3.991864886961527E-03, diy*-5.989606308601243E-03, diy*-6.441296943255076E-05); // Jupiter (mass,x,y,z,vx,vy,vz)
    planet venus(2.45E-06, -6.659025735517159E-01, 2.679456344309009E-01, 4.202093570001260E-02, diy*-7.487959643868251E-03, diy*-1.891686948899370E-02, diy*1.723630250669739E-04); // Venus (mass,x,y,z,vx,vy,vz)
    planet mercury(0.00044, -2.848623795998002E-01, -3.462144368841656E-01, -2.553868615340488E-03, diy*1.609327711354721E-02, diy*-1.647034944460669E-02, diy*-2.823041041215458E-03); // Mercury (mass,x,y,z,vx,vy,vz)
    planet neptune(0.0000515, 2.861655034326066E+01, -8.812285724589767E+00, -4.780250496232122E-01, diy* 9.027818833131188E-04, diy* 3.019033715562471E-03, diy*-8.256684941111110E-05); // Neptune (mass,x,y,z,vx,vy,vz)
    planet pluto(6.55E-09, 1.055617467780279E+01, -3.171155893980367E+01, 3.398734604938684E-01, diy*3.055434179764572E-03, diy*3.445534930791035E-04, diy*-9.194998995737056E-04); // Pluto (mass,x,y,z,vx,vy,vz)
    planet sun(1, 2.208054875983525E-03, 5.746280454272564E-03, -1.299102546588019E-04, -diy*5.245593715780954E-06, diy*5.482120330588081E-06, diy*1.232780722108486E-07);

    Solver sunearth;

    sunearth.addPlanet(earth);
    sunearth.addPlanet(mars);
    sunearth.addPlanet(saturn);
    sunearth.addPlanet(uranus);
    sunearth.addPlanet(jupiter);
    sunearth.addPlanet(venus);
    sunearth.addPlanet(mercury);
    sunearth.addPlanet(neptune);
    sunearth.addPlanet(pluto);
    sunearth.addPlanet(sun);

    vec t = linspace(0,1,10000);

    sunearth.solve(t, false);
    cout << "done";

}
