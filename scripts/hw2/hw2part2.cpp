/*----------------------------------------------------
This script was written by Lane Schultz for MSE 760
Homework 2
----------------------------------------------------*/

#include <stdio.h>                         // Print values to screen
#include <fstream>                         // Use output file

#include "lattice_fcc.cpp"                 // FCC creator
#include "force_energy_lj.cpp"             // LJ energy and accelerations
#include "reduced_units.cpp"               // Reduce units
#include "unreduced_units.cpp"             // Unreduce units
#include "velocities.cpp"                  // Random velocities

main()
{
    long double a = 5.7e-10;               // Lattice constant [m]
    long double m = 6.6e-23;               // Mass [g]
    long double T = 300.0;                 // Temperature [K]

    int n = 5;                             // Number of units cells [-]

    // Reduced units
    long double ared = reduced_units(m, 1, a);
    long double Tred = reduced_units(m, 3, T);

    // Grab the atom coordinates for FCC structure
    int dimensions;                        // Dimension of problem
    int atoms;                             // Number of atoms

    // Coordinates for FCC lattice
    long double **array = lattice_fcc(n, ared, &atoms, &dimensions);

    // The random velocity coordinates for each atom
    long double temp;  // To check temperature at end
    long double **vel = velocities(atoms, Tred, &temp);
 
    printf("Temperature: %Lf [K]", unreduced_units(m, 3, temp));
    printf("\n");
}
