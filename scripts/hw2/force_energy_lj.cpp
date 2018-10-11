/*----------------------------------------------------
This scripts calculates the cohesive energy for a 
system with the Lennard-Jones potential.
Also returns the acceleration coordinates for atoms.
----------------------------------------------------*/
#include <fstream>                         // Use output file

#include <math.h>

// Return the cohesive energy and accelerations of the system
long double force_energy_lj(
                            long double rx[],
                            long double ry[],
                            long double rz[],
                            long double ax[],
                            long double ay[],
                            long double az[],
                            long double l,
                            int         atoms,
                            int         periodic
                            )
{
    // Half lengths
    long double poshalf = l/2.0;
    long double neghalf = -poshalf;

    // The coordinates of the i atom
    long double rxi;
    long double ryi;
    long double rzi;

    // The difference between vectors
    long double drx;
    long double dry;
    long double drz;

    // The distance between atoms
    long double distance;

    // The energy between atoms
    long double u;
    long double utot = 0.0;

    // Acceleration
    long double acc;
    long double incrementax = 0.0;
    long double incrementay = 0.0;
    long double incrementaz = 0.0;

    for(int i = 0; i < atoms - 1; i++)
    {
        rxi = rx[i];
        ryi = ry[i];
        rzi = rz[i];

        for(int j = i + 1; j < atoms; j++)
        {

            switch(periodic)
            {
            case 0:                     // No periodic boundary

                drx = rxi-rx[j];
                dry = ryi-ry[j];
                drz = rzi-rz[j];

            break;

            case 1:                     // Periodic boundary

                if(rxi-rx[j] < neghalf)
                    drx = rxi-rx[j]+l;
                else if(rxi-rx[j] > poshalf)
                    drx = rxi-rx[j]-l;
                else
                    drx = rxi-rx[j];

                if(ryi-ry[j] < neghalf)
                    dry = ryi-ry[j]+l;
                else if(ryi-ry[j] > poshalf)
                    dry = ryi-ry[j]-l;
                else
                    dry = ryi-ry[j];

                if(rzi-rz[j] < neghalf)
                    drz = rzi-rz[j]+l;
                else if(rzi-rz[j] > poshalf)
                    drz = rzi-rz[j]-l;
                else
                    drz = rzi-rz[j];

            break;
            }

            distance = sqrt(pow(drx, 2.0)+pow(dry, 2.0)+pow(drz, 2.0));

            u = 1.0/pow(distance, 12.0)-1.0/pow(distance, 6.0);
            utot += u;

            acc = 1.0/pow(distance, 14.0)-0.5/pow(distance, 8.0);

            incrementax += acc*drx;
            incrementay += acc*dry;
            incrementaz += acc*drz;

            ax[i] += incrementax;
            ay[i] += incrementay;
            az[i] += incrementaz;

            ax[j] -= incrementax;
            ay[j] -= incrementay;
            az[j] -= incrementaz;

        }

        ax[i] *= 48.0;
        ay[i] *= 48.0;
        az[i] *= 48.0;

    }

    // The last atom multiplication
    ax[atoms-1] *= 48.0;
    ay[atoms-1] *= 48.0;
    az[atoms-1] *= 48.0;

    return utot*4.0;
}
