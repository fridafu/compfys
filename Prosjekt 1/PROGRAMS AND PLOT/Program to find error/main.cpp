#include <iostream>
#include <fstream>
#include "armadillo"
#include <string>
#include <iomanip>

using namespace arma;
using namespace std;




int main()
{

    int n;
    cout << "how many iterations? n = ?" << endl;
    cin >> n;


    double h = 1./(n+1);


    vec v(n);
    // open file
    ifstream inputFile("/Users/stianbilek/Documents/C++ Scripts/p1error/n" + to_string(n) + ".txt");
    string linebuffer;
    int i = 0;
    while ( getline(inputFile, linebuffer) ) {
        v(i) = stof(linebuffer);
        i++;
    }

    inputFile.close();




    double error;
    vec x = linspace<vec>(0,1,n);
    double maxval = 0;
    for (int i=0; i < n; i++)
    {


        error = abs((v[i] - (1 - (1 - exp(-10))*(x[i]+h) - exp(-10*(x[i]+h))))/v[i]);

        if(abs(error) > maxval)
        {
            maxval = abs(error);
        }


    }



    cout << scientific << setprecision(10) << maxval << endl;

    return 0;
}

/*
 *10
 *2.3511225065e+00
 * */
/*
100
2.0395723426e+00
 * */
/*
 1000
2.0039921790e+00
 * */
/*
 * 10000000
 2.0000004987e+00
 * */


