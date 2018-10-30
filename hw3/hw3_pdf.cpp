/*----------------------------------------------------
This script was written by Lane Schultz for MSE 760
Homework 2
----------------------------------------------------*/

#include <fstream>                         // Use output file

#include "lattice_fcc.cpp"                 // FCC coordinates

#include "reduced_units.cpp"               // Reduce units
#include "unreduced_units.cpp"             // Unreduce units

#include "mcmove.cpp"                      // Monte Carlo

main()
{
    long double a = 5.256e-10;             // Lattice constant [m]
    long double m = 6.6e-26;               // Mass [-/atom]
    long double sigma = 3.4e-10;           // Length [m]
    long double epsilon = 0.0104;          // Energy [eV]

    long double k = 8.6173303e-5;          // Boltzmann constant [eV/K]

    long double T = 240.0;                 // Temperature [K]

    int n = 7;                             // Number of units cells
    int atoms = n*n*n*4;                   // Number of atoms
    long double l = n*a;                   // Side length of box

    // Reduced units
    long double ared = reduced_units(m, epsilon, sigma, 1, a);
    long double lred = reduced_units(m, epsilon, sigma, 1, l);
    long double Tred = reduced_units(m, epsilon, sigma, 3, T);
	long double rhored = 0.84;             // Reduced density

    // Coordinates
    long double rx[atoms];
    long double ry[atoms];
    long double rz[atoms];

    // Coordinates for FCC lattice (here because of the needed PDF)
    lattice_fcc(n, ared, rx, ry, rz);

	// Start Monte Carlo
	mcmove(atoms, l, rx, ry, rz);
}
