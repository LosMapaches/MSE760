/*----------------------------------------------------
This scripts calculates the cohesive energy for a 
system with the Lennard-Jones potential.
----------------------------------------------------*/

#include <math.h>
#include "reduced_units.cpp"
#include "unreduced_units.cpp"

// Return the leonard jones energy of the system
long double energy_lj(
                      long double **array,
                      long double l,
                      int         periodic,
                      int         atoms
                      )
{
    // The coordinates of the i atom
    long double xi;
    long double yi;
    long double zi;

    // The coordinates of the j atom
    long double xj;
    long double yj;
    long double zj;

    // The difference between vectors
    long double dx;
    long double dy;
    long double dz;

    // The distance between atoms
    long double distance;

    // The energy between atoms
    long double u;
    long double utot;
    long double uout;

    utot = 0;
    for(int i = 0; i < atoms - 1; i++)
    {
        xi = array[i][0];
        yi = array[i][1];
        zi = array[i][2];

        for(int j = i + 1; j < atoms; j++)
        {
            xj = array[j][0];
            yj = array[j][1];
            zj = array[j][2];

            switch(periodic)
            {
            case 0:                     // No periodic boundary

                dx = xi-xj;
                dy = yi-yj;
                dz = zi-zj;

            break;

            case 1:                     // Periodic boundary

                if(xi-xj < -l/2)
                    dx = xi-xj+l;
                else if(xi-xj > l/2)
                    dx = xi-xj-l;
                else
                    dx = xi-xj;

                if(yi-yj < -l/2)
                    dy = yi-yj+l;
                else if(yi-yj > l/2)
                    dy = yi-yj-l;
                else
                    dy = yi-yj;

                if(zi-zj < -l/2)
                    dz = zi-zj+l;
                else if(zi-zj > l/2)
                    dz = zi-zj-l;
                else
                    dz = zi-zj;

            break;
            }

            distance = sqrt(pow(dx, 2)+pow(dy, 2)+pow(dz, 2));

            u = 1.0/pow(distance, 12)-1.0/pow(distance, 6);

            utot += u;
        }
    }

    uout = utot*4.0;

    return uout;
}
