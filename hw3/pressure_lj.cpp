/*----------------------------------------------------
This scripts calculates the cohesive energy for a 
system with the Lennard-Jones potential.
Also returns the acceleration coordinates for atoms.
----------------------------------------------------*/

#include <math.h>

// Return the cohesive energy and accelerations of the system
void pressure_lj(
                 long double rx[],
                 long double ry[],
                 long double rz[],
                 long double ax[],
                 long double ay[],
                 long double az[],
                 long double l,
                 long double T,
                 long double rho,
                 int         atoms,
                 int         periodic,
                 long double &pressure
                 )
{
    pressure = 0.0;

    // Half lengths
    long double poshalf = l/2.0;
    long double neghalf = -poshalf;

    // The difference between vectors
    long double drx;
    long double dry;
    long double drz;
    long double dax;
    long double day;
    long double daz;

    // The distance between atoms
    long double distance;

    // Acceleration
    long double acc;
    long double incrementax;
    long double incrementay;
    long double incrementaz;

    long double vir = 0.0;

    // Clear acceleration values
    for(int i = 0; i < atoms; i++)
    {
        ax[i] = 0.0;
        ay[i] = 0.0;
        az[i] = 0.0;
    }

    for(int i = 0; i < atoms - 1; i++)
    {
        for(int j = i + 1; j < atoms; j++)
        {
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

            acc = 48.0*(1.0/pow(distance, 14.0)-0.5/pow(distance, 8.0));

            // Limit the number of computations
            incrementax = acc*drx;
            incrementay = acc*dry;
            incrementaz = acc*drz;

            ax[i] += incrementax;
            ay[i] += incrementay;
            az[i] += incrementaz;

            ax[j] -= incrementax;
            ay[j] -= incrementay;
            az[j] -= incrementaz;

            dax = ax[i]-ax[j];
            day = ay[i]-ax[j];
            daz = az[i]-az[j];

            vir += drx*dax+dry*day+drz*daz;
        }
    }
    vir /= 3.0;
    long double volume = pow(l, 3.0); // Because NVT
    pressure = rho*T+vir/volume;
}