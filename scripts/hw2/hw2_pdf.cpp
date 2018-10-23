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

#include "pdf.cpp"                         // Pair Distribution Function

main()
{
    long double a = 5.7e-10;               // Lattice constant [m]
    long double m = 6.6e-26;               // Mass [kg]
    long double sigma = 3.4e-10;           // Length [m]
    long double epsilon = 0.0104;          // Energy [eV]
    long double rho = 1.449;               // Density [kg/m^3]

    long double k = 8.6173303e-5;          // Boltzmann constant [eV/K]

    long double T = 500.0;                 // Temperature [K]

    int n = 5;                             // Number of units cells [-]
    int atoms = n*n*n*4;                   // Number of atoms
    long double l = n*a;                   // Side length of box

    long double tred = 0.005;              // Time step
    int steps = 50;                        // The number of steps

    int bins = 15;                         // The number of bisn for RDF
    int periodic = 1;                      // Turn on PBC

    // Reduced units
    long double ared = reduced_units(m, epsilon, sigma, 1, a);
    long double lred = reduced_units(m, epsilon, sigma, 1, l);
    long double Tred = reduced_units(m, epsilon, sigma, 3, T);
    long double rhored = reduced_units(m, epsilon, sigma, 6, rho);

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

    // RDF values
    long double gr[bins];
    long double dist[bins];

    // Coordinates for FCC lattice (here because of the needed PDF)
    lattice_fcc(n, ared, rx, ry, rz);

    // Calculate the pair distribution fuction at the beginning
    pdf(rx, ry, rz, lred, rhored, atoms, periodic, bins, gr, dist);

    // Export the PDF for the inital configuration
    std::ofstream radialstart;
    radialstart.open("./pdf/radialstart");

    radialstart << "g(r) dist[A]";
    radialstart << "\n";

    for(int i = 0; i < bins; i++)
    {
        radialstart << gr[i] << " ";
        radialstart << unreduced_units(m, epsilon, sigma, 1, dist[i]) << "\n";
    }
    radialstart.close();

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
             steps,
             periodic
             );

    // Calculate the pair distribution fuction at the end
    pdf(rx, ry, rz, lred, rhored, atoms, periodic, bins, gr, dist);

    // Export the PDF for the end configuration
    std::ofstream radialend;
    radialend.open("./pdf/radialend");

    radialend << "g(r) dist[A]";
    radialend << "\n";

    for(int i = 0; i < bins; i++)
    {
        radialend << gr[i] << " ";
        radialend << unreduced_units(m, epsilon, sigma, 1, dist[i]) << "\n";
    }
    radialend.close();

    // Final Temperature
    printf("Final Temperature: %Lf [K]", unreduced_units(m, epsilon, sigma, 3, kinetic[steps]));

    // Export the energies of the system with respect to time
    std::ofstream energies;
    energies.open("./energies/energylarge");

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
