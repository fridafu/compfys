#include "ising.h"

Ising::Ising(double coupling, int l, double temp)
{
    arma_rng::set_seed_random();
    T = temp;
    beta = 1/(T*1.38064852e-23);
    stepcount = 1;

    expE = 0;
    expE2 = 0;
    expM = 0;
    expM2 = 0;
    expabsM = 0;
    L = l;
    J = coupling;
    state = rand_state();
    E = energy();
    M = magnetization();


}
void Ising::set_state(mat S)
{
    state = S;
}
mat Ising::flip_rand_spin(mat S)
{

    fr = rand() % L;
    fc = rand() % L;
    S(fr, fc) = -1*S(fr, fc);
    return S;
}
mat Ising::rand_state()
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
double Ising::magnetization()
{
    M = accu(state);
    return M;
}

double Ising::energy()
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
        E += state(L-1,i)*state(L-1,i+1);
        E += state(i,L-1)*state(i,0) + state(L-1,i)*state(0,i);
    }
    E += state(L-1,L-1)*state(L-1,0) + state(L-1,L-1)*state(0,L-1);

    E = -J*E;
    return E;
}
void Ising::step_exp_vals()
{
    for (int i = 0; i < L*L; i++)
    {
        fr = rand() % L;
        fc = rand() % L;
        up = fr + 1;
        down = fr - 1;
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
                left = L - 1;
            }
        }


        dE = 2*J*state(fr,fc)*(state(up,fc) + state(down,fc) + state(fr,right) + state(fr,left));
        dM = 2*(-state(fr,fc));

        if (dE <= 0 || double(rand())/RAND_MAX < exp(-beta*dE) )
        {
            state(fr,fc) = -state(fr,fc);
            E +=  dE;
            M += dM;
        }
    }

    E2 = E*E;
    M2 = M*M;
    absM = abs(M);

    stepcount += 1;

}
double Ising::heat_capacity()
{
    if (stepcount == 0)
    {
        cout << "Expectation values not calculated" << endl;
    }
    double k = 1.38064852e-23;
    return (E2 - E*E)/(k*T*T);
}
double Ising::magnetic_susceptibility()
{
    if (stepcount == 0)
    {
        cout << "Expectation values not calculated" << endl;
    }
    double k = 1.38064852e-23;
    return (M2 - M*M)/(k*T);
}
vec Ising::exp_vals(int steps)
{

    for (long int i = 0; i < steps; i++)
    {
        step_exp_vals();
        expE += E;
        expE2 += E*E;
        expM += M;
        expM2 += M*M;
        expabsM += absM;

    }
    E = expE/stepcount;
    E2 = expE2/stepcount;
    M = expM/stepcount;
    M2 = expM2/stepcount;
    absM = expabsM/stepcount;
    return {E, E2, M, absM, M2};
}

