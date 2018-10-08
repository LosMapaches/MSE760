/*----------------------------------------------------
This script was written by Lane Schultz for MSE 760
Homework 2
----------------------------------------------------*/
#include <stdio.h>

#include <fstream>                         // Use output file
#include "lattice_fcc.cpp"                 // FCC creator
#include "force_energy_lj.cpp"             // LJ energy

main()
{
    std::ofstream coordinates;
    coordinates.open("coordinates.txt");   // Cohesive energy vs size

    long double a = 5.7e-10;               // Lattice constant [m]
    long double m = 6.6e-23;               // Mass [g]

    int n = 2;                             // Number of units cells [-]

    // Reduce the lattice constant
    long double ared = reduced_units(m, 1, a);
    long double l = n*ared;                // Side length [m]

    // Grab the atom coordinates for FCC structure
    int dimensions;                        // Dimension of problem
    int atoms;                             // Number of atoms

    // Cohesive energy
    long double energy;

    // Set a loop to spit out cohesive energy for multiple system sizes
    coordinates << "x y z \n";

    long double **array = lattice_fcc(n, ared, &atoms, &dimensions);

    long double **acc = force_energy_lj(array, l, 1, atoms, &energy);

    for(int i=0; i<atoms; i++)
    {
        printf("%Lf", acc[i][0]);
        printf(" ");
        printf("%Lf", acc[i][1]);
        printf(" ");
        printf("%Lf", acc[i][2]);
        printf(" ");
        printf("\n");

    }

    long double ucoh = unreduced_units(m, 2, energy)/atoms;
    printf("%Lf", ucoh);

    coordinates << ucoh << "\n";
    coordinates.close();
}
