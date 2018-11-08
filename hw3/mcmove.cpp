/*----------------------------------------------------
This performs Monete Carlo simulation
----------------------------------------------------*/

#include <cstdlib>      // srand, rand
#include "energy.cpp"   // Calculate energy

void mcmove(
            int         atoms,
	    long double l,
	    long double epsilon,
	    long double T,
	    long double delta,
	    long double rx[],
	    long double ry[],
	    long double rz[],
	    long double &energyout,
	    int         &accept
	    )

{
    // Constants
    long double e = 2.718281828;
    long double beta = 1.0/(T*epsilon);

    int index = rand() % atoms;  // Index of random atom
    long double energy1 = 0.0;
    long double energy2 = 0.0;

    // Trial data
    long double trialrx[atoms];
    long double trialry[atoms];
    long double trialrz[atoms];

    // Calculate the energy of the atom
    energy(rx, ry, rz, l, atoms, index, 1, energy1);

    // Random Displacement
    long double random1 = (long double)rand()/(long double)(RAND_MAX);
    long double random2 = (long double)rand()/(long double)(RAND_MAX);
    long double random3 = (long double)rand()/(long double)(RAND_MAX);

    trialrx[index] = (long double) rx[index]+(random1-0.5)*delta;
    trialry[index] = (long double) ry[index]+(random2-0.5)*delta;
    trialrz[index] = (long double) rz[index]+(random3-0.5)*delta;

    // Assign trial
    for(int i = 0; i < atoms; i++)
    {
        // Skip the filled index
        if(i == index)
            i += 1;

        trialrx[i] = rx[i];
        trialry[i] = ry[i];
        trialrz[i] = rz[i];
    }
    // Calculate the energy of the atom
    energy(trialrx, trialry, trialrz, l, atoms, index, 1, energy2);

    long double randomcriterion = (long double)rand()/(long double)(RAND_MAX);

    // Acceptance criterion
    if(randomcriterion < pow(e, -beta*(energy2-energy1)))
    {
        rx[index] = trialrx[index];
        ry[index] = trialry[index];
        rz[index] = trialrz[index];

        energyout = energy2;
        accept = 1;
    }
    else
    {
        energyout = energy1;
        accept = 0;
    }
}
