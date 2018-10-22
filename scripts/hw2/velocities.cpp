/*----------------------------------------------------------------------------
This script randomly assigns velocities to atoms.
It is assumed that all atoms have the same velocity magnitude.
----------------------------------------------------------------------------*/

#include <cstdlib>                            // Random number generator

void velocities(
                long double vx[],
                long double vy[],
                long double vz[],
                int atoms,
                long double temperature
                )
{
    long double tempper = temperature*3.0;  // The temperature for each atom

    // Random number on unit sphere
    long double s2;
    long double zeta0;
    long double zeta1;

    // The maginitude for calulating a unit vector
    long double mag;

    int i = 0;
    while(i < atoms)
    {
        zeta0 = (long double) 2.0*rand()/RAND_MAX-1.0;
        zeta1 = (long double) 2.0*rand()/RAND_MAX-1.0;
        s2 = pow(zeta0, 2.0)+pow(zeta1, 2.0);
        if(s2 < 1)
        {
            vx[i] = 2*pow(1-s2, 0.5)*zeta0;
            vy[i] = 2*pow(1-s2, 0.5)*zeta1;
            vz[i] = 1-2*s2;

            mag = pow(pow(vx[i], 2.0)+pow(vy[i], 2.0)+pow(vz[i], 2.0), 0.5);

            vx[i] *= tempper/mag;
            vy[i] *= tempper/mag;
            vz[i] *= tempper/mag;

            i++;  // Increment the while loop
        }
    }
}
