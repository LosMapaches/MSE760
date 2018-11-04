/*----------------------------------------------------
This scripts calculates the cohesive energy for a 
system with the Lennard-Jones potential.
Also returns the acceleration coordinates for atoms.
----------------------------------------------------*/

#include <math.h>

// Return the cohesive energy and accelerations of the system
void energy(
            long double rx[],
            long double ry[],
            long double rz[],
            long double l,
            int         atoms,
            int         i,
            int         periodic,
            long double &cohesive
            )
{
    // Half lengths
    long double poshalf = l/2.0;
    long double neghalf = -poshalf;

    // The difference between vectors
    long double drx;
    long double dry;
    long double drz;

    // The distance between atoms
    long double distance;

    // The energy between atoms
    long double u;
    cohesive = 0.0;

    for(int j = 0; j < atoms; j++)
    {
        // Skip the starting atom
	if(j == i)
	    j += 1;

        switch(periodic)
        {
        case 0:                     // No periodic boundary
            drx = rx[i]-rx[j];
            dry = ry[i]-ry[j];
            drz = rz[i]-rz[j];

        break;

        case 1:                     // Periodic boundary
            if(rx[i]-rx[j] < neghalf)
                drx = rx[i]-rx[j]+l;
            else if(rx[i]-rx[j] > poshalf)
                drx = rx[i]-rx[j]-l;
            else
                drx = rx[i]-rx[j];

            if(ry[i]-ry[j] < neghalf)
                dry = ry[i]-ry[j]+l;
            else if(ry[i]-ry[j] > poshalf)
                dry = ry[i]-ry[j]-l;
            else
                dry = ry[i]-ry[j];

            if(rz[i]-rz[j] < neghalf)
                drz = rz[i]-rz[j]+l;
            else if(rz[i]-rz[j] > poshalf)
                drz = rz[i]-rz[j]-l;
            else
                drz = rz[i]-rz[j];

        break;
        }
        distance = sqrt(pow(drx, 2.0)+pow(dry, 2.0)+pow(drz, 2.0));
        u = 1.0/pow(distance, 12.0)-1.0/pow(distance, 6.0);
        cohesive += u;
    }
    cohesive *= 4.0;
}
