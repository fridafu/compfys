#include <iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include<time.h>
#include<algorithm>
#include<vector>
#include "jacobi.h"
#include <string>

using namespace std;
using namespace arma;

//Creating a tri-diagonal matrix with diagonal elements d+V, neighbouring elements to the diagonal e, and the rest of the off diagonal elements = 0
mat makeA(double rho_min, double rho_max, int n, bool interact, double wr) {
    // Step size
    double h = (rho_max - rho_min)/(n);
    // Init potential
    double V = 0;
    // Making a matrix filled with zeros
    mat A = mat(n,n); A.zeros();
    double d = 2./(h*h);
    double e = -1./(h*h);
    int i=0;
    //iterating over values of rho
    for (i = 0; i < n; i++) {
        double rho = (i+1)*h;
        // Coulomb potential (interaction) or not?
        if (interact) {
            V = wr*wr*rho*rho + 1./(rho);
            //cout << "YEY" <<endl;
        } else {
            V = rho*rho;
            //cout << "NOOO" << endl;
        }

       // Assigning elements on the diagonal with the potential
       A(i,i) = d + V;
       // Editing the side-diagonal elements
       if (i < n-1){
           A(i,i+1) = e;
           A(i+1,i) = e;
        }
    }
    return A;

}
// performs jacobi algorithm
// to find eigenvalues/vectors
int jacobi(int n, int interact, double conv, double wr, mat& a, mat& v) {
    cout.precision(5);
    double aip=0, aiq=0, vpi=0, vqi=0;
    double tau=0, t=0, s=0, c=0;//tan(theta), sin(theta), cos(theta)
    int count=1;                //count of iterations
    int count_old=count-10;     //keep track of every 10th iteration
    int p=n-1, q=n-2;           //off diag all same value to start
                                //pick last as first maximum
    clock_t start, end;

    if(n<=10){
        cout<<"Before diagonalization"<<endl;
        print_vals(a,v,n,conv);
        cout<<endl;
    }

    double app=a(p,p);
    double aqq=a(q,q);
    double apq=a(p,q);

    start=clock();

    while(abs(apq)>conv){
        if(count>1){
            apq=0;
            find_max(a,p,q,apq,n);
        }

        //calculate sin(theta) and cos(theta)
        aqq=a(q,q);
        app=a(p,p);
        tau=(aqq-app)/(2*apq);
        if(tau>0)
            t=1/(tau+sqrt(1+tau*tau));
        else
            t=-1/(-tau+sqrt(1+tau*tau));
        c=1/sqrt(1+t*t);
        s=c*t;

        //calculate new matrix elements and vectors
        for(int i=0;i<n;i++){
            if(i!=p && i!=q){
                aip=a(i,p);
                aiq=a(i,q);
                a(i,p)=aip*c-aiq*s;
                a(p,i)=aip*c-aiq*s;
                a(i,q)=aiq*c+aip*s;
                a(q,i)=aiq*c+aip*s;
            }
            //vpi=v(p,i);
            //vqi=v(q,i);
            vpi=v(i,p);
            vqi=v(i,q);
            //v(p,i)=c*vpi-s*vqi;
           // v(q,i)=c*vqi+s*vpi;
            v(i,p)=c*vpi-s*vqi;
            v(i,q)=c*vqi+s*vpi;
        }
        a(p,p)=app*c*c-2*apq*c*s+aqq*s*s;
        a(q,q)=app*s*s+2*apq*c*s+aqq*c*c;
        a(p,q)=0;
        a(q,p)=0;

        count++;
    }

    end=clock();

    if(n<=10){
        cout<<"After diagonalization"<<endl;
        print_vals(a,v,n,conv);
        cout<<endl;
    }

    cout<<"Diagonalization took "<<count<<" iterations"<<endl;
    cout<<scientific<<"CPU time (sec) : "<<((double)end-(double)start)/CLOCKS_PER_SEC<<endl;

    return 0;
}

//get first three eigenvectors
mat get_eigenvecs(mat a, mat v, int n, bool write_to_file){
    vector<double>eigenvals=get_eigenvals(a,n);
    mat vecs(3,n);
    for(int i=0;i<3;i++){
        for(int j=0;j<n;j++){
            if(a(j,j)==eigenvals[i]){
                for(int k=0;k<n;k++){
                      vecs(i,k)=v(k,j);
                }
             }
         }
    }

    if (write_to_file) {
        // writing eigenvectors to file
        string filename;
        cout << "Give me a filname: ";
        cin >> filename;
        ofstream myfile;
        myfile.open("/home/hannahcb/compfys/Prosjekt2/" +  filename +".txt");
        myfile << vecs << endl;
        myfile.close();
    }

    return vecs;
}

//get eigenvalues in order
vector<double> get_eigenvals(mat a,int n){
    vector<double>eigen;
    for(int i=0;i<n;i++){
        eigen.push_back(a(i,i));
    }
    sort (eigen.begin(), eigen.begin()+n);
    return eigen;
}

//find maximum non-diag matrix elements
void find_max(mat a,int& p,int& q,double& apq,int n){
    apq = 0;
    for (int i=0;i<n;i++){
         for (int j=0;j<n;j++){
            if(i!=j && abs(a(i,j))>=abs(apq)){
                apq=a(i,j);
                p=i;
                q=j;
            }
         }
    }
}

//print matrix and eigenvectors
void print_vals(mat A, mat v,int n,double conv){
    cout<<"A: ";
    for (int i=0;i<n;i++){
        if(i>0){
            cout<<"   ";
         }
        for (int j=0;j<n;j++){
            if(abs(A(i,j))>conv)
                cout<<fixed<<A(i,j)<<" ";
            else cout<<"0.000 ";
        }
        cout<<endl;
    }
    for (int i=0;i<n;i++){
        cout<<"v"<<i<<": ";
        for (int j=0;j<n;j++){
            if(abs(v(j,i))>conv)
                cout<<fixed<<v(j,i)<<" ";
            else cout<<"0.000 ";
        }
        cout<<endl;
    }
}

void test() {

    try{

        int q, r;
        int n = 100;
        double maxoffdiag;

        mat B = { {1, 3, 0, 0, 0},
                  {2, 4, 6, 0, 0},
                  {0, 3, 7, 2, 0},
                  {0, 0, 7, 9, 4},
                  {0, 0, 0, 1, 3}
                };

       find_max(B, q, r, maxoffdiag, 5);

       /*Checks if find_max() returns correct value and indexes
        */
       if (maxoffdiag != 7 || q != 3 || r != 2) {
           throw 99;

       }

       mat A = makeA(0, 5, n, false, 0); //Makes matrix with known eigenvalues
       mat v = eye(n,n);

       //Disables cout during test.

       streambuf* orig_buf = cout.rdbuf();
       cout.rdbuf(NULL);
       int jac = jacobi(n, true, 1e-8, 0, A, v);
       cout.rdbuf(orig_buf);

       vec ei = get_eigenvals(A, n);

       /*
        * Throw exception if the relative error in eigenvalues is more than 0.1%
        */
       for (int u = 0; u < n; u++) {
           if ( (1 - ei(u)/( 2*(2*u + 3./2.) ) )*100  > 0.1) {
               throw 'a';
               break;
           }
       }


       mat first_three_vectors = get_eigenvecs(A, v, n, false);

       double eps = 1e-8;

       /* Calculates the inner product of
        * the eigenvectors and throws exception
        * it they are not orthogonal.
       */
       for (int i = 0; i < 3; i++) {

           for (int j = 0; j < 3; j++) {

               double inner_product = 0;

               for (int u = 0; u < n; u++) {
                    inner_product += first_three_vectors(i,u)*first_three_vectors(j,u);
               }
               if (i != j && fabs( inner_product ) > eps) {
                   throw 1.5;
                   break;
               }
               if (i == j && fabs(fabs( inner_product ) - 1 ) > eps) {
                   throw 1.5;
                   break;
               }
           }
       }


       }
       catch(int x) {
           cout << "ERROR: maxdiag() did not return the correct matrix element!"<<endl;
       }
       catch(char a) {
           cout <<"ERROR: jacobi() did not return the correct eigenvalues!"<<endl;
       }
       catch(double b) {
           cout <<"ERROR: get_eigenvecs() did not return orthogonal eigenvectors!"<<endl;
       }
}


int main(){
    cout << "Unit test initiated!" << endl;
    test(); //run unit tests
    cout << "Unit test complete!" << endl;

    //Boundary "close enough" to 0
    double eps = 1e-8;
    double rho_min = 0;
    double rho_max;
    //"frequency" wr reflects coulomb potential
    double wr;
    double n;
    cout << "Gimme an n: " ;
    cin >> n;
    cout << "Gimme a rho_max: (5 is good) ";
    cin >> rho_max;
    cout << "Gimme an wr: " ;
    cin >> wr;
    //Make a matrix with the given parameters
    mat Amat = makeA(rho_min, rho_max, n, true, wr);
   //Create identity matrix v
    mat v = mat(n,n); v.zeros();
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            if(i==j){
                v(i,j)=1;
            }
            else{
                v(i,j)=0;
            }
        }
    }
    //Run the jacobi algorithm for our matrix and parameters
    int jac = jacobi(n, true, eps, wr, Amat, v);
    //Finding the first three eigenvectors of the sorted vectors.
    mat first_three_vectors = get_eigenvecs(Amat, v, n, true);
    cout << endl;
    first_three_vectors.print();
    // finding the first three eigenvalues
    vec lam = get_eigenvals(Amat,n);
    cout << lam(0) << endl << lam(1) << endl << lam(2) << endl;
    //Taking the time it takes to run the Armadillo solver
    clock_t start2, end2;
    start2 = clock();
    //Solving the eigenvalue problem with a built in Armadillo solver
    vec ADsolver = eig_sym(Amat);
    end2 = clock();
    cout<<scientific<<"Armadillo CPU time (sec) : "<<((double)end2-(double)start2)/CLOCKS_PER_SEC<<endl;
    cout << "Armadillo says = " << endl << ADsolver(0) << endl << ADsolver(1) << endl << ADsolver(2) << endl;
    return 0;
}
