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
            long double rxatom,
            long double ryatom,
            long double rzatom,
            int         index,
            long double l,
            int         atoms,
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
	    i++;

        switch(periodic)
        {
        case 0:                     // No periodic boundary
            drx = rx[i]-rxatom;
            dry = ry[i]-ryatom;
            drz = rz[i]-rzatom;

        break;

        case 1:                     // Periodic boundary
            if(rx[i]-rxatom < neghalf)
                drx = rx[i]-rxatom+l;
            else if(rx[i]-rxatom > poshalf)
                drx = rx[i]-rxatom-l;
            else
                drx = rx[i]-rxatom;

            if(ry[i]-ryatom < neghalf)
                dry = ry[i]-ryatom+l;
            else if(ry[i]-ryatom > poshalf)
                dry = ry[i]-ryatom-l;
            else
                dry = ry[i]-ryatom;

            if(rz[i]-rzatom < neghalf)
                drz = rz[i]-rzatom+l;
            else if(rz[i]-rzatom > poshalf)
                drz = rz[i]-rzatom-l;
            else
                drz = rz[i]-rzatom;

        break;
        }
        distance = sqrt(pow(drx, 2.0)+pow(dry, 2.0)+pow(drz, 2.0));
        energy += 1.0/pow(distance, 12.0)-1.0/pow(distance, 6.0);
    }
    energy *= 4.0;
}
