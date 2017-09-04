#include <iostream>

#include <fstream>

#include "armadillo"


using namespace arma;
using namespace std;

vec gaussian_elim(vec a, vec b, vec c, vec btilde) {

    int n = size(b)(0);

    vec v(n);
    vec dtilde(n);
    vec ytilde(n);

    dtilde(0) = b(0);
    ytilde(0) = btilde(0);

    for (int i = 1; i < n; i++) {
        dtilde(i) = b(i) - c(i-1)*c(i-1)/dtilde(i-1);
        ytilde(i) = btilde(i) - c(i-1)*ytilde(i-1)/dtilde(i-1);
    }

    v(n-1) = ytilde(n-1)/dtilde(n-1);

    for (int j = n-2; j >= 1; j--) {
        v(j) = ( ytilde(j) - c(j)*v(j+1)/dtilde(j+1))/dtilde(j);
    }

    return v;
}

vec gaussian_special(int n, vec btilde, float a, float b, float c) {
    vec v(n);
    float ytilde = btilde(0);
    float ytilde1 = 0;  //ytilde1 er ytilde(i+1). ytilde er ytilde(i)
    v(n-1) = ytilde(n-1)/(-(n)/(n-1.));
    for (int i = n-2; i >= 0; i--) {
        ytilde1 = btilde(i+1) - c*ytilde/(-(i + 1.)/i);
        v(i) = (ytilde - c*v(i+1)/(-(i+2.)/(i+1.)));
        ytilde = ytilde1;

    }

    return v;

}



int main()
{
    int n;
    cin >> n;

    float h = 1./(n+1);
    /*
    vec a = -1*ones<vec>(n-1);
    vec b = 2*ones<vec>(n);
    vec c = -1*ones<vec>(n-1);
    */
    float a = -1;
    float b = 2;
    float c = -1;
    vec x = linspace<vec>(0, 1, n);
    vec f = h*h*100*exp(-10*x);

    vec data = gaussian_special(n, f, a, b, c);

    cout << data << endl;

    return 0;
}
