#include <iostream>

#include <fstream>

#include "armadillo"


using namespace arma;
using namespace std;

vec gaussian_elim(vec a, vec b, vec c, vec d) {

    int n = size(b)(0);

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



int main()
{
    int n;
    cin >> n;

    float h = 1./(n+1);

    vec a = -1*ones<vec>(n-1);
    vec b = 2*ones<vec>(n);
    vec c = -1*ones<vec>(n-1);

    vec x = linspace<vec>(0, 1, n);
    vec f = h*h*100*exp(-10*x);

    //vec data = gaussian_elim(a, b, c, f);
    vec data = gaussian_special(n, f);

    cout << data << endl;

    return 0;
}
