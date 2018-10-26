/*----------------------------------------------------
This script reduces the units of the problem
----------------------------------------------------*/

#include <math.h>                           // square root

// Return the reduced unit
long double reduced_units(
                          long double sigma,
                          long double epsilon,
                          long double m,
                          int   choice,
                          long double i 
                          )
{
    /* Constants for conversion */
    long double k = 1.38e-23;              // Boltzmann constant [J/K]

    long double ri;                        // The reduced output

    switch(choice)
    {
        case 1: ri = i/sigma; // Convert distance
            break;
        case 2: ri = i/epsilon; // Convert energy
            break;
        case 3: ri = i/(epsilon/k); // Convert temperature
            break;
        case 4: ri = i/sqrt((m*pow(sigma, 2))/epsilon); // Convert time
            break;
        case 5: ri = i/(epsilon/sigma); // Convert force
            break;
    }

    return ri;
}
