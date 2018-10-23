/*----------------------------------------------------
This script reduces the units of the problem
----------------------------------------------------*/

#include <math.h>                           // square root

// Return the reduced unit
long double reduced_units(
                          long double m,
                          long double epsilon,
                          long double sigma,
                          int         choice,
                          long double i
                          )
{
    long double k = 8.6173303e-5;           // Boltzmann constant [eV/K]

    switch(choice)
    {
        case 1: i /= sigma; // Convert distance
            break;
        case 2: i /= epsilon; // Convert energy
            break;
        case 3: i /= epsilon/k; // Convert temperature
            break;
        case 4: i /= sqrt((m*pow(sigma, 2.0))/epsilon); // Convert time
            break;
        case 5: i /= epsilon/sigma; // Convert force
            break;
        case 6: i *= pow(sigma, 3.0); // Convert deinsty
            break;
    }
    return i;
}
