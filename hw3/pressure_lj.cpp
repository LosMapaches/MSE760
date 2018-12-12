/*----------------------------------------------------
This scripts calculates the pressure for a 
system with the Lennard-Jones potential.
Also returns the acceleration coordinates for atoms.
----------------------------------------------------*/

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

    // Vector magnitudes
    long double magr;
    long double maga;

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

            dax = ax[i]-ax[j];
            day = ay[i]-ax[j];
            daz = az[i]-az[j];

            magr = sqrt(drx*drx+dry*dry+drz*drz);
            maga = sqrt(dax*dax+day*day+daz*daz);

            pressure += magr*maga;
        }
    }
    pressure /= 3.0;
    pressure /= pow(l, 3.0);
    pressure += rho*T;
}
