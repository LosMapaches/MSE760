/*----------------------------------------------------
This script was written by Lane Schultz for MSE 760
Homework 2
----------------------------------------------------*/

#include <fstream>                         // Use output file

#include "lattice_fcc.cpp"                 // FCC creator
#include "force_energy_lj.cpp"             // LJ energy and accelerations
#include "reduced_units.cpp"               // Reduce units
#include "unreduced_units.cpp"             // Unreduce units

main()
{
    std::ofstream coordinates;
    coordinates.open("coordinates.txt");   // Cohesive energy vs size

    long double a = 5.7e-10;               // Lattice constant [m]
    long double m = 6.6e-23;               // Mass [g]

    int n = 5;                             // Number of units cells [-]

    // Reduced units
    long double ared = reduced_units(m, 1, a);
    long double l = n*ared;                // Side length

    // Grab the atom coordinates for FCC structure
    int dimensions;                        // Dimension of problem
    int atoms;                             // Number of atoms

    // Coordinates for FCC lattice
    long double **array = lattice_fcc(n, ared, &atoms, &dimensions);

    // The acceleration coordinates for each atom
    long double energy;  // Where energy is stored in reduced units
    long double **acc = force_energy_lj(array, l, 1, atoms, &energy);

    coordinates << "atom x[N] y[N] z[N]";
    coordinates << "\n";
    for(int i = 0; i < atoms; i++)
    {
        coordinates << i << " ";
        coordinates << unreduced_units(m, 5, acc[i][0])*1.602e-19 << " ";
        coordinates << unreduced_units(m, 5, acc[i][1])*1.602e-19 << " ";
        coordinates << unreduced_units(m, 5, acc[i][2])*1.602e-19 << " ";
        coordinates << "\n";

    }

    coordinates.close();
}
