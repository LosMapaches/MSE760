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
#include "temperature.cpp"                 // Calculate temperature

#include "simulation.cpp"                  // Start molecular dynamics

main()
{
    long double a = 5.7e-10;               // Lattice constant [m]
    long double m = 6.6e-26;               // Mass [kg]
    long double sigma = 3.4e-10;           // Length [m]
    long double epsilon = 0.0104;          // Energy [eV]
    long double k = 8.6173303e-5;          // Boltzmann constant [eV/K]

    long double T = 300.0;                 // Temperature [K]

    int n = 5;                             // Number of units cells [-]
    int atoms = n*n*n*4;                   // Number of atoms
    long double l = n*a;                   // Side length of box

    /*
    long double tred = 0.001;              // Time step
    int steps = 10000;                     // The number of steps
    */

    /*
    long double tred = 0.005;              // Time step
    int steps = 2000;                      // The number of steps
    */

    long double tred = 0.02;               // Time step
    int steps = 500;                       // The number of steps

    // Reduced units
    long double ared = reduced_units(m, epsilon, sigma, 1, a);
    long double lred = reduced_units(m, epsilon, sigma, 1, l);
    long double Tred = reduced_units(m, epsilon, sigma, 3, T);

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

    // The energies for each step calulated
    long double cohesive[steps];
    long double kinetic[steps];
    long double total[steps];

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

    // Export the energies of the system with respect to time
    std::ofstream energies;
    energies.open("./energies/energies0p02");

    energies << "time[s] cohesive[eV/atom] kinetic[eV/atom] total[eV/atom]";
    energies << "\n";

    for(int i = 0; i <= steps; i++)
    {
        energies << unreduced_units(m, epsilon, sigma, 4, tred*i) << " ";

        cohesive[i] = unreduced_units(m, epsilon, sigma, 2, cohesive[i]);
        cohesive[i] /= atoms;
        energies << cohesive[i] << " ";

        kinetic[i] = unreduced_units(m, epsilon, sigma, 3, kinetic[i]);
        kinetic[i] *= 3.0/2.0*k;  // Normalized by atoms
        energies << kinetic[i] << " ";

        energies << cohesive[i]+kinetic[i] << "\n";
    }
    energies.close();
}
