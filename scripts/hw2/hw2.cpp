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
#include "simulation.cpp"                  // Simulate

main()
{
    std::ofstream coordinates;
    coordinates.open("coordinates.txt");   // Cohesive energy vs size

    long double a = 5.7e-10;               // Lattice constant [m]
    long double m = 6.6e-23;               // Mass [g]
    long double T = 300.0;                 // Temperature [K]
    long double t = 2.2e-11;               // Time [s]
    long double timestep = 0.001e-12;      // Time step [s/step]

    int n = 5;                             // Number of units cells [-]

    // Reduced units
    long double ared = reduced_units(m, 1, a);
    long double tred = reduced_units(m, 4, t);
    long double timestepred = reduced_units(m, 4, timestep);
    long double Tred = reduced_units(m, 3, T);
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

    // The random velocity coordinates for each atom
    long double temp;  // To check temperature at end
    long double **vel = velocities(atoms, Tred, &temp);
 
    printf("Temperature: %Lf [K]", unreduced_units(m, 3, temp));
    printf("\n");
    
    // The cohesive energy of the system
    long double ucoh = unreduced_units(m, 2, energy)/atoms;
    printf("Cohesive Energy: %Lf [eV]", ucoh);
    printf("\n");

    int steps = 10; // The number of run steps
    long double **r = simulate(atoms, steps, l, Tred, timestepred, array);

    coordinates.close();
}
