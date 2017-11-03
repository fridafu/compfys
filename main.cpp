#include <iostream>
#include <armadillo>
#include <random>
#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace arma;
using namespace std;

string outfilename = "p4.txt";
double k = 1.38064852e-23;
double initial_temp = 1.0;
double final_temp = 3.0;
double temp_step = 1e-6;
double J = 1.; //coupling constant
vec T = linspace(initial_temp,final_temp,temp_step);
vec beta = 1./(k*T);
int n_spins = 2.;
double E = 0;
double M = 0;
int mcs = 1e6; //monte carlo cycles
vec w(17);
vec average(5);
ofstream ofile;

//2X2-case to compare with

vec Z = 2*exp(-8*J*beta) + 2*exp(8*J*beta) + 12;

vec E_exp = (J/2.)*(16*exp(8*J*beta) - 16*exp(-8*J*beta));

// inline function for PeriodicBoundary boundary conditions
   inline int PeriodicBoundary(int i, int limit, int add) {
        return (i+limit+add) % (limit);
        }

// Functions
void initialize(int, double, mat, double, double);
void Metropolis(int, long&, mat, double, double, vec); // prints to file the results of the calculations
void output(int, int, double, vec);

// function to initialise energy, spin matrix and magnetization
void initialize(int n_spins, mat spin_matrix,  double E, double M)
{
    // setup spin matrix and initial magnetization
  for(int x =0; x < n_spins; x++) {
    for (int y= 0; y < n_spins; y++){
      spin_matrix(x,y) = 1.0; // spin orientation for the ground state
      M +=  (double) spin_matrix(x,y);
    }
  }

  // setup initial energy
  for(int x =0; x < n_spins; x++) {
    for (int y= 0; y < n_spins; y++){
      E -=  (double) spin_matrix(x,y)*(spin_matrix(PeriodicBoundary(x,n_spins,-1),y) +
     spin_matrix(x,PeriodicBoundary(y,n_spins,-1)));
    }
  }
}


int main(){
        long idum;
        mat spin_matrix = mat(n_spins, n_spins);
        ofile.open(outfilename);
        idum = -1; // random starting point

        for ( double temp = initial_temp; temp <= final_temp; temp+=temp_step){
            // initialise energy and magnetization
            E = M = 0.;
            // setup array for possible energy changes
            for( int de =-8; de <= 8; de++) w(de+8) = 0;
            for( int de =-8; de <= 8; de+=4) w(de+8) = exp(-de/temp); // initialise array for expectation values
            for( int i = 0; i < 5; i++) average(i) = 0.;
            initialize(n_spins, spin_matrix, E, M);

            // start Monte Carlo computation
            for (int cycles = 1; cycles <= mcs; cycles++){
                Metropolis(n_spins, idum, spin_matrix, E, M, w);
        // update expectation values
        average(0) += E; average(1) += E*E; average(2) += M; average(3) += M*M; average(4) += fabs(M);
        // print results
        output(n_spins, mcs, temp, average);
        }
    }
    ofile.close(); // close output file return 0;
}


void Metropolis(int n_spins, long& idum, mat spin_matrix, double E, double M, vec w) {
    // loop over all spins
    for(int y =0; y < n_spins; y++){
        for (int x= 0; x < n_spins; x++){
            // Find random position
            int ix = (int) (rand()*(double)n_spins);
            int iy = (int) (rand()*(double)n_spins);
            int deltaE = 2*spin_matrix(ix,iy)*
                    (spin_matrix(ix,PeriodicBoundary(iy,n_spins,-1))+
                     spin_matrix(PeriodicBoundary(ix,n_spins,-1),iy) +
                     spin_matrix(ix,PeriodicBoundary(iy,n_spins,1)) +
                     spin_matrix(PeriodicBoundary(ix,n_spins,1),iy));

            // Metropolis test
            if (rand() <= w(deltaE+8) ) {
                spin_matrix(iy,ix) *= -1; // flip one spin and accept new spin config

            // update energy and magnetization
            M += (double) 2*spin_matrix(iy,ix);
            E += (double) deltaE;
    }
}
}
}



//Wanna print to file?
void output(int n_spins, int mcs, double temperature, vec average) {
    double norm = 1/((double) (mcs)); // divided by total number of cycles
    double Eaverage = average[0]*norm;
    double E2average = average[1]*norm;
    double Maverage = average[2]*norm;
    double M2average = average[3]*norm;
    double Mabsaverage = average[4]*norm;
    // all expectation values are per spin, divide by 1/n_spins/n_spins
    double Evariance = (E2average- Eaverage*Eaverage)/n_spins/n_spins;
    double Mvariance = (M2average - Maverage*Maverage)/n_spins/n_spins;
    double M2variance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << temperature;
    ofile << setw(15) << setprecision(8) << Eaverage/n_spins/n_spins;
    ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature; // ofile << setw(15) << setprecision(8) << Maverage/n_spins/n_spins;
    ofile << setw(15) << setprecision(8) << M2variance/temperature;
    ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins << endl;
} // end output function
