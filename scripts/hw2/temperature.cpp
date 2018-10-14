/*----------------------------------------------
Calculate the temperature with velocities
----------------------------------------------*/

#include <math.h>                           // square root

long double temperature(
                        int atoms,
                        long double vx[],
                        long double vy[],
                        long double vz[]
                        )

{

    // A Check on temperature
    long double temp = 0.0;

    for(int i = 0; i < atoms; i++)
    {
        temp += pow(pow(vx[i], 2.0)+pow(vy[i], 2.0)+pow(vz[i], 2.0), 0.5);        
    }

    temp /= 3.0*atoms;  // The temperature

    return temp;
}
