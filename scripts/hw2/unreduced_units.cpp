/*----------------------------------------------------
This script unreduces the units of the problem
----------------------------------------------------*/

#include <math.h>                           // square root

// Return standard units 
long double unreduced_units(
                            long double m,
                            int         choice,
                            long double ri 
                            )
{
    /* Constants for conversion */
    long double k = 1.38e-23;              // Boltzmann constant [J/K]
    long double sigma = 3.4e-10;           // Length [m]
    long double epsilon = 0.0104;          // Energy [eV]

    long double i;                         // The reduced output

    switch(choice)
    {
        case 1: i = ri*sigma; // Convert distance
            break;
        case 2: i = ri*epsilon; // Convert energy
            break;
        case 3: i = ri*(epsilon/k); // Convert temperature
            break;
        case 4: i = ri*sqrt((m*pow(sigma, 2))/epsilon); // Convert time
            break;
        case 5: i = ri*(epsilon/sigma); // Convert force
            break;
    }

    return i;
}
