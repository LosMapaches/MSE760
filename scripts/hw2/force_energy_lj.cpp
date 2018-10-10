/*----------------------------------------------------
This scripts calculates the cohesive energy for a 
system with the Lennard-Jones potential.
Also returns the acceleration coordinates for atoms.
----------------------------------------------------*/

#include <math.h>

// Return the leonard jones energy of the system
long double **force_energy_lj(
                              long double **array,
                              long double l,
                              int         periodic,
                              int         atoms,
                              long double *energy
                              )
{
    // The coordinates of the i atom
    long double xi;
    long double yi;
    long double zi;

    // Half lengths
    long double poshalf = l/2.0;
    long double neghalf = -poshalf;

    // The difference between vectors
    long double dx;
    long double dy;
    long double dz;

    // The distance between atoms
    long double distance;

    // The energy between atoms
    long double u;
    long double utot;

    // Acceleration
    long double a;
    long double ax = 0.0;
    long double ay = 0.0;
    long double az = 0.0;

    // Matrix to hold acceleration values
    long double **acc = new long double *[atoms];
    for(size_t i = 0; i < atoms; i++)
        acc[i] = new long double [3];

    utot = 0.0;
    for(int i = 0; i < atoms - 1; i++)
    {
        xi = array[i][0];
        yi = array[i][1];
        zi = array[i][2];

        for(int j = i + 1; j < atoms; j++)
        {

            switch(periodic)
            {
            case 0:                     // No periodic boundary

                dx = xi-array[j][0];
                dy = yi-array[j][1];
                dz = zi-array[j][2];

            break;

            case 1:                     // Periodic boundary

                if(xi-array[j][0] < neghalf)
                    dx = xi-array[j][0]+l;
                else if(xi-array[j][0] > poshalf)
                    dx = xi-array[j][0]-l;
                else
                    dx = xi-array[j][0];

                if(yi-array[j][1] < neghalf)
                    dy = yi-array[j][1]+l;
                else if(yi-array[j][1] > poshalf)
                    dy = yi-array[j][1]-l;
                else
                    dy = yi-array[j][1];

                if(zi-array[j][2] < neghalf)
                    dz = zi-array[j][2]+l;
                else if(zi-array[j][2] > poshalf)
                    dz = zi-array[j][2]-l;
                else
                    dz = zi-array[j][2];

            break;
            }

            distance = pow(pow(dx, 2.0)+pow(dy, 2.0)+pow(dz, 2.0), 0.5);

            u = 1.0/pow(distance, 12.0)-1.0/pow(distance, 6.0);
            utot +=u;

            a = 1.0/pow(distance, 14.0)-0.5/pow(distance, 8.0);

            ax += a*dx;
            ay += a*dy;
            az += a*dz;

            acc[i][0] += ax;
            acc[i][1] += ay;
            acc[i][2] += az;

            acc[j][0] -= ax;
            acc[j][1] -= ay;
            acc[j][2] -= az;

        }

        acc[i][0] *= 48.0;
        acc[i][1] *= 48.0;
        acc[i][2] *= 48.0;

    }

    // The last atom multiplication
    acc[atoms-1][0] *= 48.0;
    acc[atoms-1][1] *= 48.0;
    acc[atoms-1][3] *= 48.0;

    *energy = utot*4.0;

    return acc;

    // Delete pointer
    delete energy;

    // Delete the matrix
    for(size_t i = 0; i < atoms; i++)
        delete acc[i];
    delete acc;

}
