/*----------------------------------------------------
This performs Monete Carlo simulation
----------------------------------------------------*/

#include "particle_energy.cpp"             // Atom Energy

void mcmove(
            int         atoms,
	    long double l,
	    long double T,
	    long double delta,
	    long double rx[],
	    long double ry[],
	    long double rz[],
            long double &cohesive,
            int         periodic,
	    long double &energyout,
	    int         &accept
	    )

{
    int index = rand() % atoms;  // Index of random atom
    long double energy1;
    long double energy2;

    // Calculate the energy of the atom
    particle_energy(rx, ry, rz, rx[index], ry[index], rz[index], index, l, atoms, periodic, energy1);

    // Random Numbers [0, 1]
    long double random1 = (long double)rand()/(long double)(RAND_MAX);
    long double random2 = (long double)rand()/(long double)(RAND_MAX);
    long double random3 = (long double)rand()/(long double)(RAND_MAX);
    long double randomcriterion = (long double)rand()/(long double)(RAND_MAX);

    // Trial move
    long double trialrx = rx[index]+(random1-0.5)*delta;
    long double trialry = ry[index]+(random2-0.5)*delta;
    long double trialrz = rz[index]+(random3-0.5)*delta;

    // Calculate the energy of the atom
    particle_energy(rx, ry, rz, trialrx, trialry, trialrz, index, l, atoms, periodic, energy2);

    // Acceptance criterion
    if(randomcriterion < (long double) exp(-(energy2-energy1)/T))
    {
        rx[index] = trialrx;
        ry[index] = trialry;
        rz[index] = trialrz;

        energyout = energy2;
        accept = 1;
    }

    else
    {
        energyout = energy1;
        accept = 0;
    }
    cohesive += energyout-energy1;
    energyout = cohesive/atoms;
}
