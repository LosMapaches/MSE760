/*----------------------------------------------------
This script was written by Lane Schultz for MSE 760
Homework 2
----------------------------------------------------*/

#include <fstream>                         // Use output file

#include "lattice_fcc.cpp"                 // FCC coordinates

#include "reduced_units.cpp"               // Reduce units
#include "unreduced_units.cpp"             // Unreduce unitsi

#include "force_energy_lj.cpp"             // LJ energy and accelerations
#include "velocities.cpp"                  // Randomize velocities

#include "simulation.cpp"                  // Start molecular dynamics

main()
{
    long double a = 5.7e-10;                  // Lattice constant [m]
    long double m = 6.6e-23;                  // Mass [g]
    long double sigma = 3.4e-10;              // Length [m]
    long double epsilon = 0.0104;             // Energy [eV]

    long double T = 300.0;                    // Temperature [K]

    int n = 5;                                // Number of units cells [-]
    int atoms = n*n*n*4;                      // Number of atoms
    long double l = n*a;                      // Side length of box

    long double t = 0.001e-12;                // Time step [s/step]
    int steps = 22000;                        // The number of steps dt = 0.001 [ps]

    // Reduced units
    long double ared = reduced_units(m, epsilon, sigma, 1, a);
    long double lred = reduced_units(m, epsilon, sigma, 1, l);
    long double Tred = reduced_units(m, epsilon, sigma, 3, T);
    long double tred = reduced_units(m, epsilon, sigma, 4, t);

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

    // Start simulation
    steps = 22000;  // dt = 0.001

    // The energies for each step calulated
    long double cohesive[steps];
    long double kinetic[steps];
    long double total[steps];

    // Temperatures
    long double temp[steps];
    simulate(
             atoms,
             n,
             ared,
             Tred,
             cohesive,
             kinetic,
             rx,
             ry,
             rz,
             vx,
             vy,
             vz,
             ax,
             ay,
             az,
             lred,
             tred,
             steps
             );

    steps = 4400;   // dt = 0.005
    steps = 1100;   // dt = 0.02
}
