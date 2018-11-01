/*----------------------------------------------------
This performs Monete Calro simulation
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

	srand(time(NULL));  // Generate random seed

	int index = rand() % atoms;  // Index of random atom
	long double cohesive1;
	long double cohesive2;

	// Trial data
	long double trialrx[atoms];
	long double trialry[atoms];
	long double trialrz[atoms];

	// Calculate the energy of the atom
	energy(rx, ry, rz, l, atoms, index, 1, cohesive1);

	// Random Displacement
	trialrx[index] = rx[index]+(rand()%2-0.5)*delta;
    trialry[index] = ry[index]+(rand()%2-0.5)*delta;
    trialrz[index] = rz[index]+(rand()%2-0.5)*delta;

	// Calcualte trial energy
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
	energy(trialrx, trialry, trialrz, l, atoms, index, 1, cohesive2);

	// Acceptance criterion
	if(rand()%2 < pow(e, -beta*(cohesive2-cohesive1)))
	{
		for(int i = 0; i < atoms; i++)
		{
			rx[i] = trialrx[i];
            ry[i] = trialry[i];
            rz[i] = trialrz[i];
		}
		energyout = cohesive2;
		accept = 1;
	}
	else
		energyout = cohesive1;
		accept = 0;
}
