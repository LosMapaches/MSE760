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
    r_coordinates.open("r_coordinates.txt");  // Position Coordinates

    std::ofstream a_coordinates;
    a_coordinates.open("a_coordinates.txt");  // Acceleration Coordinates


    long double a = 5.7e-10;                  // Lattice constant [m]
    long double m = 6.6e-23;                  // Mass [g]
    long double k = 1.38e-23;                 // Boltzmann constant [J/K]
    long double sigma = 3.4e-10;              // Length [m]
    long double epsilon = 0.0104;             // Energy [eV]

    int n = 5;                                // Number of units cells [-]

    // Reduced units
    reduced_units(m, epsilon, sigma, 1, a);
    long double l = n*a;                      // Side length

    // Grab the atom coordinates for FCC structure
    int atoms = n*n*n*4;                      // Number of atoms

    // Coordinates
    long double rx[atoms];
    long double ry[atoms];
    long double rz[atoms];

    // Velocities
    long double vx[atoms];
    long double vy[atoms];
    long double vz[atoms];

    // Accelerations
    long double ax[atoms];
    long double ay[atoms];
    long double az[atoms];

    // Coordinates for FCC lattice
    lattice_fcc(n, a, rx, ry, rz);

    r_coordinates << "atom x[m] y[m] z[m]";
    r_coordinates << "\n";
    for(int i = 0; i < atoms; i++)
    {
        unreduced_units(m, epsilon, sigma, 1, rx[i]);
        unreduced_units(m, epsilon, sigma, 1, ry[i]);
        unreduced_units(m, epsilon, sigma, 1, rz[i]);

        r_coordinates << i << " ";
        r_coordinates << rx[i] << " ";
        r_coordinates << ry[i] << " ";
        r_coordinates << rz[i] << " ";
        r_coordinates << "\n";

    }

    // The acceleration coordinates for each atom
    long double energy;  // Where energy is stored in reduced units
    force_energy_lj(rx, ry, rz, ax, ay, az, energy, l, 1, atoms);

    a_coordinates << "atom x[N] y[N] z[N]";
    a_coordinates << "\n";
    for(int i = 0; i < atoms; i++)
    {
        unreduced_units(m, epsilon, sigma, 5, ax[i]);
        unreduced_units(m, epsilon, sigma, 5, ay[i]);
        unreduced_units(m, epsilon, sigma, 5, az[i]);

        a_coordinates << i << " ";
        a_coordinates << ax[i]*1.602e-19 << " ";
        a_coordinates << ay[i]*1.602e-19 << " ";
        a_coordinates << az[i]*1.602e-19 << " ";
        a_coordinates << "\n";

    }

    r_coordinates.close();
    a_coordinates.close();
}
