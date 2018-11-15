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
            int         index,
            int         periodic,
            long double &energy
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
    energy = 0.0;

    for(int i = 0; i < atoms; i++)
    {
        // Skip the starting atom
	if(index == i)
	    i += 1;

        switch(periodic)
        {
        case 0:                     // No periodic boundary
            drx = rx[i]-rx[index];
            dry = ry[i]-ry[index];
            drz = rz[i]-rz[index];

        break;

        case 1:                     // Periodic boundary
            if(rx[i]-rx[index] < neghalf)
                drx = rx[i]-rx[index]+l;
            else if(rx[i]-rx[index] > poshalf)
                drx = rx[i]-rx[index]-l;
            else
                drx = rx[i]-rx[index];

            if(ry[i]-ry[index] < neghalf)
                dry = ry[i]-ry[index]+l;
            else if(ry[i]-ry[index] > poshalf)
                dry = ry[i]-ry[index]-l;
            else
                dry = ry[i]-ry[index];

            if(rz[i]-rz[index] < neghalf)
                drz = rz[i]-rz[index]+l;
            else if(rz[i]-rz[index] > poshalf)
                drz = rz[i]-rz[index]-l;
            else
                drz = rz[i]-rz[index];

        break;
        }
        distance = sqrt(pow(drx, 2.0)+pow(dry, 2.0)+pow(drz, 2.0));
        energy += 1.0/pow(distance, 12.0)-1.0/pow(distance, 6.0);
    }
    energy *= 4.0;
}
