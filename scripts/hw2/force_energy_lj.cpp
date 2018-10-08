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

    // Acceleration
    long double a;
    long double ax = 0;
    long double ay = 0;
    long double az = 0;

    // Matrix to hold acceleration values
    long double **acc = new long double *[atoms];
    for(size_t i = 0; i < atoms; i++)
        acc[i] = new long double [3];

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
            utot +=u;

            a = 1.0/pow(distance, 14)-0.5/pow(distance, 8);

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

    // Delete the matrix
    for(size_t i = 0; i < atoms; i++)
        delete acc[i];
    delete acc;

}
