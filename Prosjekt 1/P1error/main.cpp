#include <iostream>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

int const n = 10000000;
int const np = n + 1;
double h = 1./(n+1);
vec x = linspace<vec>(0,1,n);

//double h[6] = {1, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7};
//double exact = 1./3;
//double f(double x)
//{
//    double ex = 100*exp(10*x));
//    return ex;
//}
//double u(double x)
//    double uex = 1 - (1-100*exp(-10))*x - exp(-10*x);
//    return uex;

int main()
{
    //ifstream v;
    //v.open("tall.txt");
    //cout << v << endl;

    vec v(n);
    // open file
    ifstream inputFile("/home/frida/Documents/FYS4150/P1error/n10000000.txt");
    string linebuffer;
    int i = 0;
    while ( getline(inputFile, linebuffer) ) {
        v(i) = stof(linebuffer);
        i++;
    }

    cout << v << endl;
    inputFile.close();
    //double df2[n];
    //double maxval;
    //double error2[n];
    ofstream outFile("/home/frida/Documents/FYS4150/P1error/error10000000.txt");
    double u[n];
    double error[n];
    for (int i=0; i < n; i++)
    {
        //double v[n];
        //df2[i] = ((f(x+(h[i])) - f(x)))/h[i];
        //cout << "Derivative of arctan(" << x << ") = " << df2[i] << ", with h = " << h[i] << endl;
        //error2[i] = (abs(df2[i] - exact)/exact);
        //cout << "Error = " << error2[i] << endl;

        u[i] = 1 - (1-100*exp(-10))*(x[i]+h) - exp(-10*(x[i]+h));
        error[i] = abs((v[i]-u[i])/v[i]);
        outFile << error[i] << endl;



        cout << "Error = " << error[i] << endl;
    }

    outFile.close();
    double maxval = 0;
    for (int i=0; i<n; i++)
    {
        if(abs(error[i]) > maxval)
        {
            maxval = abs(error[i]);
        }
    }
    cout << "Maximum error: " << maxval << endl;
    return 0;
}
