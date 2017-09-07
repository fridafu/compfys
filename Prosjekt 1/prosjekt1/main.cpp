#include <iostream>
#include <fstream>
#include "armadillo"
#include <ctime>
#include <string>

using namespace arma;
using namespace std;

vec gaussian_elim(vec a, vec b, vec c, vec d) {

    int n = size(b,0);

    vec v(n);
    vec btilde(n);
    vec dtilde(n);
    vec dtt(n);

    btilde(0) = b(0);
    dtilde(0) = d(0);

    for (int i = 1; i < n; i++) {
        btilde(i) = b(i) - a(i-1)*c(i-1)/btilde(i-1);
        dtilde(i) = d(i) - a(i-1)*dtilde(i-1)/btilde(i-1);
    }

    v(n-1) = dtilde(n-1)/btilde(n-1);
    dtt(n-1) = dtilde(n-1);

    for (int j = n-2; j >= 0; j--) {
        dtt(j) = dtilde(j) - c(j)*dtt(j+1)/btilde(j+1);
        v(j) = dtt(j)/btilde(j);
    }

    return v;
}

vec gaussian_special(int n, vec d) {

    vec v(n);
    vec dtt(n);
    vec dtilde(n);

    dtilde(0) = d(0);

    for (int i = 1; i<n; i++) {
        dtilde(i) = d(i) + dtilde(i-1)/((i+1.)/i);

    }

    v(n-1) = dtilde(n-1)/((n+1.)/n);
    dtt(n-1) = dtilde(n-1);

    for (int j = n-2; j >= 0; j--) {
        dtt(j) = dtilde(j) + dtt(j+1)/((j+3.)/(j+2));
        v(j) = dtt(j)/((j+2.)/(j+1));
    }

    return v;

}

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

    int n;

    cin >> n;

    float h = 1./(n+1);

    vec a = -1*ones<vec>(n-1);
    vec b = 2*ones<vec>(n);
    vec c = -1*ones<vec>(n-1);
    vec x = linspace<vec>(0, 1, n);
    vec f = h*h*100*exp(-10*(x+h));
    clock_t t;

    t = clock();
    vec data = gaussian_elim(a, b, c, f);
    t = clock() - t;

    cout << "gaussian_elim t = " << double(t)/CLOCKS_PER_SEC << " s" << endl;

    t = clock();
    vec data1 = gaussian_special(n, f);
    t = clock() - t;

    cout << "gaussian_special t = " << double(t)/CLOCKS_PER_SEC << " s" << endl;

    t = clock();
    vec data2 = solvethri(n, -1, 2, -1, f);
    t = clock() - t;

    cout << "solvethri t = " << double(t)/CLOCKS_PER_SEC  << " s " << "total time of 8 CPUs" << endl; //divided by 8

    ofstream myfile;

    myfile.open("/Users/stianbilek/Documents/Python Scripts/n" + to_string(n) +".txt");

    for (int i = 0; i < n; i++) {
        myfile << data(i) << endl;
    }

    myfile.close();


    return 0;
}

/*
 c)
 *1000000                                       //input
gaussian_elim t = 0.157942 s
gaussian_special t = 0.063577 s
*/

/*
 e
 *10000                                 //10^4 worked, but not 10^5
gaussian_elim t = 0.001734 s
gaussian_special t = 0.000794 s
solvethri t = 54.9899 s total time of 8 CPUs               //divide by 8 to get time we had to wait
*/


