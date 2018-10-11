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
    long double a = 5.7e-10;               // Lattice constant [m]
    long double m = 6.6e-23;               // Mass [g]
    long double T = 300.0;                 // Temperature [K]
    long double timestep = 0.001e-12;      // Time step [s/step]
    int steps;

    int n = 2;                             // Number of units cells [-]

    // Reduced units
    long double ared = reduced_units(m, 1, a);
    long double timestepred = reduced_units(m, 4, timestep);
    long double Tred = reduced_units(m, 3, T);
    long double l = n*ared;                // Side length

    // Run for 2.2e-11 seconds
    steps = 22000;  // The number of run steps when dt = 0.001 [ps]
    long double **r = simulate(n, ared, l, Tred, timestepred, steps);

    steps = 4400;  // The number of run steps when dt = 0.005 [ps]
    steps = 1100;  // The number of run steps when dt = 0.02 [ps]
}
