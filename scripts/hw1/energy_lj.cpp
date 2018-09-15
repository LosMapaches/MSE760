/*----------------------------------------------------
This script creates an fcc lattice in 3D
----------------------------------------------------*/

#include <math.h>
#include "reduced_units.cpp"
#include "unreduced_units.cpp"

// Return the leonard jones energy of the system
long double energy_lj(
                      long double **array,
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
    long double U;
    long double Utot = 0;
    long double Uout;

    for(int i = 0; i < atoms-1; i++)
    {
        for(int j = i + 1; j < atoms; j++)
        {
            xi = reduced_units(sigma, epsilon, m, 1, array[i][0]);
            yi = reduced_units(sigma, epsilon, m, 1, array[i][1]);
            zi = reduced_units(sigma, epsilon, m, 1, array[i][2]);

            xj = reduced_units(sigma, epsilon, m, 1, array[j][0]);
            yj = reduced_units(sigma, epsilon, m, 1, array[j][1]);
            zj = reduced_units(sigma, epsilon, m, 1, array[j][2]);

            dx = xi-xj;
            dy = xi-xj;
            dz = xi-xj;

            distance = sqrt(pow(dx, 2)+pow(dy, 2)+pow(dz, 2));

            U = 4*(1/pow(distance, 12)-1/pow(distance, 6));

            if(isnan(U) == 1)
                Utot += 0;
            else
                Utot += U;

        }
    }

    Uout = unreduced_units(sigma, epsilon, m, 2, Utot);

    return Uout;
}
