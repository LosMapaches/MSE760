/*----------------------------------------------------------------------------
This script randomly assigns velocities to atoms.
It is assumed that all atoms have the same velocity magnitude.
----------------------------------------------------------------------------*/

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

    // Random number on unit sphere
    long double s;
    long double zeta0;
    long double zeta1;
    long double xi0 [atoms];
    long double xi1 [atoms];
    long double xi2 [atoms];

    int i = 0;
    while(i < atoms)
    {
        zeta0 = (long double) 2.0*rand()/RAND_MAX-1.0;
        zeta1 = (long double) 2.0*rand()/RAND_MAX-1.0;
        s = pow(zeta0, 2.0);
        s += pow(zeta1, 2.0);
        if(s < 1)
        {
            xi0[i] = 2*pow(1-pow(s, 2.0), 0.5)*zeta0;
            xi1[i] = 2*pow(1-pow(s, 2.0), 0.5)*zeta1;
            xi2[i] = 1-2*pow(s, 2);

            i++;
        }
    }

    // The maginitude for calulating a unit vector
    long double mag;

    // Assign random vector directions to each atom
    long double temp = 0.0;
    for(int i = 0; i < atoms; i++)
    {
        vel[i][0] = xi0[i];
        vel[i][1] = xi1[i];
        vel[i][2] = xi2[i];

        mag = pow(pow(vel[i][0], 2.0)+pow(vel[i][1], 2.0)+pow(vel[i][2], 2.0), 0.5);

        for(int j = 0; j < 3; j++)
            vel[i][j] *= tempper/mag;

        temp += pow(pow(vel[i][0], 2.0)+pow(vel[i][1], 2.0)+pow(vel[i][2], 2.0), 0.5);
    }

    temp /= 3.0*atoms;

    *tempcheck = temp;

    return vel;

    // Delete pointer
    delete tempcheck;

    // Delete the matrix
    for(size_t i = 0; i < atoms; i++)
        delete [] vel[i];
    delete [] vel;

}
