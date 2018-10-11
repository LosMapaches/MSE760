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
    std::ofstream r_coordinates;
    r_coordinates.open("r_coordinates.txt");  // Cohesive energy vs size

    long double a = 5.7e-10;                  // Lattice constant [m]
    long double m = 6.6e-23;                  // Mass [g]
    long double k = 1.38e-23;              // Boltzmann constant [J/K]
    long double sigma = 3.4e-10;           // Length [m]
    long double epsilon = 0.0104;          // Energy [eV]

    int n = 5;                                // Number of units cells [-]

    // Reduced units
    a = reduced_units(m, epsilon, sigma, 1, a);
    long double l = n*a;                   // Side length

    // Grab the atom coordinates for FCC structure
    int atoms = n*n*n*4;                      // Number of atoms

    // Coordinates for FCC lattice
    long double x[atoms];
    long double y[atoms];
    long double z[atoms];
    lattice_fcc(n, a, x, y, z);

    r_coordinates << "atom x[N] y[N] z[N]";
    r_coordinates << "\n";
    for(int i = 0; i < atoms; i++)
    {
        r_coordinates << i << " ";
        r_coordinates << unreduced_units(m, epsilon, sigma, 1, x[i]) << " ";
        r_coordinates << unreduced_units(m, epsilon, sigma, 1, y[i]) << " ";
        r_coordinates << unreduced_units(m, epsilon, sigma, 1, z[i]) << " ";
        r_coordinates << "\n";

    }

    /*
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
    */

    r_coordinates.close();
}
