/*----------------------------------------------------
This scripts calculates the cohesive energy for a 
system with the Lennard-Jones potential.
----------------------------------------------------*/

#include <stdio.h>

#include <math.h>
#include "reduced_units.cpp"
#include "unreduced_units.cpp"

// Return the leonard jones energy of the system
long double energy_lj(
                      long double **array,
                      long double l,
                      int         periodic,
                      int         atoms,
                      long double sigma,
                      long double epsilon,
                      long double m
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

    long double rl = reduced_units(sigma, epsilon, m, 1, l);

    utot = 0;
    for(int i = 0; i < atoms - 1; i++)
    {
        for(int j = i + 1; j < atoms; j++)
        {
            xi = reduced_units(sigma, epsilon, m, 1, array[i][0]);
            yi = reduced_units(sigma, epsilon, m, 1, array[i][1]);
            zi = reduced_units(sigma, epsilon, m, 1, array[i][2]);

            xj = reduced_units(sigma, epsilon, m, 1, array[j][0]);
            yj = reduced_units(sigma, epsilon, m, 1, array[j][1]);
            zj = reduced_units(sigma, epsilon, m, 1, array[j][2]);

            switch(periodic)
            {
            case 0:                     // No periodic boundary

                dx = xi-xj;
                dy = yi-yj;
                dz = zi-zj;

            break;

            case 1:                     // Periodic boundary

                if(xi-xj < -rl/2)
                    dx = xi-xj+rl;
                else if(xi-xj > rl/2)
                    dx = xi-xj-rl;
                else
                    dx = xi-xj;

                if(yi-yj < -rl/2)
                    dy = yi-yj+rl;
                else if(yi-yj > rl/2)
                    dy = yi-yj-rl;
                else
                    dy = yi-yj;

                if(zi-zj < -rl/2)
                    dz = zi-zj+rl;
                else if(zi-zj > rl/2)
                    dz = zi-zj-rl;
                else
                    dz = zi-zj;

            break;
            }

            distance = sqrt(pow(dx, 2)+pow(dy, 2)+pow(dz, 2));

            u = 4*(1/pow(distance, 12)-1/pow(distance, 6));

            utot += u;
        }
    }

    uout = unreduced_units(sigma, epsilon, m, 2, utot);

    return uout;
}
