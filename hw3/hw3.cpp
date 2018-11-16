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

    int n = 4;                             // Number of units cells
    int atoms = n*n*n*4;                   // Number of atoms
    long double l = n*a;                   // Side length of box
    int steps = 1e5;                       // Number of simulation steps
    long double beta = 1.0/(T*k);

    // Reduced units
    long double ared = reduced_units(m, epsilon, sigma, 1, a);
    long double lred = reduced_units(m, epsilon, sigma, 1, l);
    long double Tred = reduced_units(m, epsilon, sigma, 3, T);
    long double rhored = 0.84;             // Reduced density
    long double betared = reduced_units(m, epsilon, sigma, 2, 1.0/beta);
    betared = pow(betared, -1.0);

    // Coordinates
    long double rx[atoms];
    long double ry[atoms];
    long double rz[atoms];

    // Energy
    long double energyout = 0.0;

    // Atom displacement
    long double delta = 0.9*ared/lred;  // Beginning displacement criterion

    // If move is accepted
    int accept;
    int control1 = 0;  // For the summation of acceptances
    long double control2;  // For determining the percent acceptance

    long double energyoutev = 0.0;

    // Coordinates for FCC lattice
    lattice_fcc(n, ared, rx, ry, rz);

    srand(time(NULL));  // Generate random seed

    // Start Monte Carlo
    std::ofstream energies;
    energies.open("./energies.txt");

    energies << "Step AcceptanceRate Energy[eV]\n";
    int periodic = 1;  // Turn on periodic boundry conditions
    for(int i = 1; i <= steps; i++)
    {
       	mcmove(atoms, lred, betared, Tred, delta, rx, ry, rz, periodic, energyout, accept);
        energyoutev = unreduced_units(m, epsilon, sigma, 2, energyout);

	control1 += accept;
        control2 = (long double) control1/i;

        printf("Step: %i |", i);
        printf("Acceptance: %Lf |", control2);
        printf("Energy [eV]: %Lf \n", energyoutev);

        energies << i << " ";
        energies << control2 << " ";
        energies << energyoutev << "\n";
    }
}
