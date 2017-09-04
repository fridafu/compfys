#include <iostream>
#include "armadillo"
#include <fstream>

using namespace arma;
using namespace std;


vec solvethri(int x, float a, float b, float c, vec btilde) {

    int j = 0;
    mat A(x,x);
    A(0,0) = b;
    A(0,1) = c;
    A(x-1,x-2) = a;
    A(x-1,x-1) = b;
    for (int i = 1; i < x - 1; i++) {
        A(i,j) = a;
        A(i, j+1) = b;
        A(i,j+2) = c;
        j++;
    }
    return A.i()*btilde;
}


int main()
{
    float a = -1;
    float b = 2;
    float c = -1;
    int n;
    cin >> n;

    float h = 1./(n+1);

    vec x = linspace<vec>(0, 1, n);
    vec f = h*h*100*exp(-10*x);

    vec dta = solvethri(n, a, b, c, f);
    cout << dta << endl;

    /* ofstream myfile;
    myfile.open("nlik1000.txt");
    for (int i = 0; i < size(dta)(0); i++) {
        myfile << dta(i) << endl;
    }
    myfile.close();
    */

    return 0;


}
