/*----------------------------------------------------
This performs Monete Calro simulation
----------------------------------------------------*/

#include <cstdlib>      // srand, rand
#include "energy.cpp"   // Calculate energy

void mcmove(
			int         atoms,
			long double l,
			long double rx[],
			long double ry[],
			long double rz[]
			)

{
	srand(time(NULL));  // Generate random seed

	int index = rand() % atoms;  // Index of random atom
	long double cohesive;

	// Trial data
	long double trialrx[atoms];
	long double trialry[atoms];
	long double trialrz[atoms];

	// calculate the energy of the atom
	energy(rx, ry, rz, l, atoms, index, 1, cohesive);

	// Random Displacement
	long double delta = l/atoms;
	trialrx[index] = rx[index]+(rand()%2-0.5)*delta;
    trialry[index] = ry[index]+(rand()%2-0.5)*delta;
    trialrz[index] = rz[index]+(rand()%2-0.5)*delta;
}
