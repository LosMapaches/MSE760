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
    int steps;                                // The number of steps

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

    // Coordinates for FCC lattice
    lattice_fcc(n, ared, rx, ry, rz);

    // Export the position coordinates
    std::ofstream r_coordinates;
    r_coordinates.open("r_coordinates.txt");

    r_coordinates << "atom x[m] y[m] z[m]";
    r_coordinates << "\n";
    for(int i = 0; i < atoms; i++)
    {
        r_coordinates << i << " ";
        r_coordinates << unreduced_units(m, epsilon, sigma, 1, rx[i]) << " ";
        r_coordinates << unreduced_units(m, epsilon, sigma, 1, ry[i]) << " ";
        r_coordinates << unreduced_units(m, epsilon, sigma, 1, rz[i]) << " ";
        r_coordinates << "\n";
    }

    r_coordinates.close();

    // The acceleration coordinates for each atom
    long double energy_cohesive = force_energy_lj(
                                                  rx,
                                                  ry,
                                                  rz,
                                                  ax,
                                                  ay,
                                                  az,
                                                  lred,
                                                  atoms,
                                                  1
                                                  );

    // Unreduce energy
    energy_cohesive = unreduced_units(
                                      m,
                                      epsilon,
                                      sigma,
                                      2,
                                      energy_cohesive
                                      );

    printf("Cohesive Energy: %Lf \n", energy_cohesive/atoms);

    // Export the acceleration coordinates
    std::ofstream a_coordinates;
    a_coordinates.open("a_coordinates.txt");

    a_coordinates << "atom x[N] y[N] z[N]";
    a_coordinates << "\n";
    for(int i = 0; i < atoms; i++)
    {
        a_coordinates << i << " ";
        a_coordinates << unreduced_units(m, epsilon, sigma, 5, ax[i])*1.602e-19 << " ";
        a_coordinates << unreduced_units(m, epsilon, sigma, 5, ay[i])*1.602e-19 << " ";
        a_coordinates << unreduced_units(m, epsilon, sigma, 5, az[i])*1.602e-19 << " ";
        a_coordinates << "\n";
    }

    a_coordinates.close();

    long double tempcheck = velocities(vx, vy, vz, atoms, Tred);
    printf("Temperature: %Lf \n", unreduced_units(m, epsilon, sigma, 3, tempcheck));

    steps = 22000;  // dt = 0.001
    simulate(atoms, rx, ry, rz, vx, vy, vz, ax, ay, az, lred, tred, tred);

    steps = 4400;   // dt = 0.005
    steps = 1100;   // dt = 0.02
}
