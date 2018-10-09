/*----------------------------------------------------------------------------
This script randomly assigns velocities to atoms.
It is assumed that all atoms have the same velocity magnitude.
----------------------------------------------------------------------------*/
#include <stdio.h>

#include <cstdlib>                            // Random number generator

long double **velocities(int atoms, long double temperature, long double *tempcheck)
{

    long double tempper = temperature*(3.0);  // The temperature for each atom

    // Create a matrix for velocites of each atom
    long double **vel = new long double *[atoms];
    for(size_t i = 0; i < atoms; i++)
    {
        vel[i] = new long double [3];
    }

    // The maginitude for calulating a unit vector
    long double mag;

    // Assign random vector directions to each atom
    long double temp = 0.0;
    for(int i = 0; i < atoms; i++)
    {
        for(int j = 0; j < 3; j++)
            vel[i][j] = (long double) rand()/RAND_MAX;

        mag = pow(pow(vel[i][0], 2.0)+pow(vel[i][1], 2.0)+pow(vel[i][2], 2.0), 0.5);

        for(int j = 0; j < 3; j++)
            vel[i][j] *= tempper/mag;

        temp += pow(pow(vel[i][0], 2.0)+pow(vel[i][1], 2.0)+pow(vel[i][2], 2.0), 0.5);
    }

    temp /= 3.0*atoms;

    *tempcheck = temp;

    return vel;

    // Delete the matrix
    for(size_t i = 0; i < atoms; i++)
        delete vel[i];
    delete vel;

}
