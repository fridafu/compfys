#include "ising.h"

Ising::Ising(double coupling, int l, double temp)
//coupling is a constant expressing the strength of the interaction between neighboring spins
//l is the dimension of the 2D lattice (lxl)
//temp is the temperature in kelvin.
//When initializing object, a random state is generated.
{
    arma_rng::set_seed_random();
    k = 1.38064852e-23;
    T = temp;
    beta = 1./(T*k);
    stepcount = 0;
    expE = 0;
    expE2 = 0;
    expM = 0;
    expM2 = 0;
    expabsM = 0;
    L = l;
    J = coupling;
    totaccept = 0;
    state = rand_state();
    energy();
    magnetization();



}
vec Ising::get_expectation_values()
//returns expectation values {mean energy, mean energy squared, mean magnetization, mean magnetization squared, mean absolute magnetization}
{
    double EE = expE/stepcount;
    double EE2 = expE2/stepcount;
    double MM = expM/stepcount;
    double MM2 = expM2/stepcount;
    double absMM = expabsM/stepcount;
    vec exp_val(5); exp_val(0) = EE; exp_val(1) = EE2; exp_val(2) = MM; exp_val(3) = MM2; exp_val(4) = absMM;
    return exp_val;
    //return {EE,EE2,MM,MM2,absMM};
}
void Ising::reset_expectation_values()
//resets the expectation values
{
    totaccept = 0;
    stepcount = 0;
    expE = 0;
    expE2 = 0;
    expM = 0;
    expM2 = 0;
    expabsM = 0;
}
void Ising::set_state(mat S)
//takes a matrix S as input and sets it as the state.
{
    state = S;
    energy();
    magnetization();
}


mat Ising::flip_rand_spin(mat S)
//flips random spin in input matrix.
{

    fr = rand() % L;
    fc = rand() % L;
    S(fr, fc) = -1*S(fr, fc);
    return S;
}
mat Ising::rand_state()
//generates and returns a random state
{

    mat S = randu(L,L);
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            if (S(i,j) <= 0.5)
            {
                S(i,j) = -1;
            }
            else
            {
                S(i,j) = 1;
            }
        }
    }
    return S;

}
void Ising::magnetization()
//calculates the magnetization and stores it in variable M
{
    M = accu(state);
}

void Ising::energy()
//calculates the energy and stores it in variable E
{
    E = 0;
    for (int i = 0; i < L-1; i++)
    {
        for (int j = 0; j < L-1; j++)
        {
            E += state(i,j)*state(i,j+1) + state(i,j)*state(i+1,j);
        }
    }
    for (int i = 0; i < L-1; i++)
    {
        E += state(L-1,i)*state(L-1,i+1) + state(i,L-1)*state(i+1,L-1);
        E += state(i,L-1)*state(i,0) + state(L-1,i)*state(0,i);
    }
    E += state(L-1,L-1)*state(L-1,0) + state(L-1,L-1)*state(0,L-1);

    E = -J*E;
}

double Ising::get_magnetization()
{
    return M;
}
double Ising::get_energy()
//returns the energy of the system
{
    return E;
}

void Ising::step_metropolis()
//performs one step of the metropolis algorithm
//stores the number of times this function is used in variable stepcount
//stores the number of accepted configurations in variable totaccept
{
    for (int i = 0; i < L*L; i++)
    {

        fr = rand() % L;
        fc = rand() % L;

        up = fr - 1;
        down = fr + 1;
        right = fc + 1;
        left = fc - 1;


        if (fr == 0)
        {
            up = L - 1;

            if (L == 2)
            {
                down = 1;
            }

        }
        else if (fr == L - 1)
        {
            down = 0;
            if (L == 2)
            {
                up = 0;
            }
        }
        if (fc == 0)
        {
            left = L - 1;
            if (L == 2)
            {
                right = L - 1;
            }
        }
        else if (fc == L - 1)
        {
            right = 0;

            if (L == 2)
            {
                left = 0;
            }
        }



        dE = 2*J*state(fr,fc)*(state(up,fc) + state(down,fc) + state(fr,right) + state(fr,left));
        dM = 2*(-state(fr,fc));

        r = double(rand())/RAND_MAX;


        if ( (dE <= 0) ||  (r <= exp(-beta*dE)) )
        {
            state(fr,fc) = -state(fr,fc);
            totaccept += 1;
            E +=  dE;
            M += dM;
        }


    }

    stepcount += 1;

}

int Ising::get_configurations()
//returns the total accepted configurations.
{
    return totaccept;
}

double Ising::heat_capacity()
//returns the heat capacity
{
    if (stepcount == 0)
    {
        cout << "Expectation values not calculated" << endl;
    }

    return (expE2/stepcount - (expE/stepcount)*(expE/stepcount))/(T*T*k);
}
double Ising::magnetic_susceptibility()
//returns the magnetic susceptibility
{
    if (stepcount == 0)
    {
        cout << "Expectation values not calculated" << endl;
    }

    return (expM2/stepcount - expabsM*expabsM/(stepcount*stepcount))/(T*k);
}
void Ising::exp_vals(int steps)
//calculates expectation values with steps steps.
{

    for (long int i = 0; i < steps; i++)
    {
        step_metropolis();
        expE += E;
        expE2 += E*E;
        expM += M;
        expM2 += M*M;
        expabsM += abs(M);
    }
}
double Ising::energy_probability(int initial_steps, int steps, int En)
//This function will do 'inital_steps' steps with the metropolis algorithm, check the energy of the system
//and then calculate the probability of the this energy occuring.
{
    exp_vals(initial_steps);
    reset_expectation_values();
    long int counter = 0;
    for (long int i = 0; i < steps; i++)
    {
        exp_vals(1);
        if (En == E)
        {
            counter++;
        }

    }
    return double(counter)/steps;
}

