/*----------------------------------------------------
This script was written by Lane Schultz for MSE 760
Homework 2
----------------------------------------------------*/

#include <fstream>                         // Use output file
#include <math.h>                          // Powers

#include "lattice_fcc.cpp"                 // FCC coordinates
#include "reduced_units.cpp"               // Reduce units
#include "unreduced_units.cpp"             // Unreduce units
#include "force_energy_lj.cpp"             // Accelerations
#include "mcmove.cpp"                      // Monte Carlo

main()
{
    long double m = 6.6e-26;               // Mass [atom^-1]
    long double sigma = 3.4e-10;           // Length [m]
    long double epsilon = 0.0104;          // Energy [eV]
    long double k = 8.6173303e-5;          // Boltzmann constant [eV/K]

    int n = 5;                             // Number of units cells
    int atoms = n*n*n*4;                   // Number of atoms
    int steps = 1e6;                       // Number of simulation steps

    // Reduced units
    long double T = 2.0;
    long double rho = 0.84;             // Reduced density
    long double l = pow(n/rho, 1.0/3.0);
    long double a = l/n;

    int periodic = 1;  // Turn on periodic boundry conditions

    // Coordinates
    long double rx[atoms];
    long double ry[atoms];
    long double rz[atoms];

    // Accelerations
    long double ax[atoms];
    long double ay[atoms];
    long double az[atoms];

    // Energy
    long double energyout = 0.0;
    long double cohesive = 0.0;

    // Pressure
    long double P = 0.0;
    long double pressuresum = 0.0;

    // Atom displacement
    long double delta = 0.2;  // Beginning displacement criterion

    // If move is accepted
    int accept = 0;
    int control1 = 0;  // For the summation of acceptances
    long double control2 = 0.0;  // For determining the percent acceptance
    long double dummy = 0.0;  // Dummy variable

    // Coordinates for FCC lattice
    lattice_fcc(n, atoms, a, rx, ry, rz);

    // Calculate the cohesive energy
    force_energy_lj(
                    rx,
                    ry,
                    rz,
                    ax,
                    ay,
                    az,
                    l,
                    T,
                    rho,
                    atoms,
                    periodic,
                    cohesive,
                    P
                    );

    printf("%Lf", unreduced_units(m, epsilon, sigma, 2, cohesive)/atoms);

    // Monte Carlo Parameters
    int frequency = 100;  // The data acquisition rate
    int count = 0;  // The number of times data is taken
    int stepstart = 400000;  // Approximate number of steps after settling

    // Start Monte Carlo
    srand(time(NULL));  // Generate random seed

    std::ofstream energies;
    energies.open("./energies.txt");

    std::ofstream pressures;
    pressures.open("./pressures.txt");

    energies << "Step AcceptanceRate Energy[eV]\n";
    pressures << "Pressure \n";

    // for(int i = 1; i <= steps; i++)
    for(int i = 1; i <= 1; i++)
    {
       	mcmove(
               atoms,
               l,
               T,
               delta,
               rx,
               ry,
               rz,
               periodic,
               cohesive,
               energyout,
               accept
               );

        energyout = unreduced_units(m, epsilon, sigma, 2, energyout);

	control1 += accept;
        control2 = (long double) control1/i;

        printf("Step: %i |", i);
        printf("Acceptance: %Lf |", control2);
        printf("Energy [eV/atom]: %Lf \n", energyout);

        energies << i << " ";
        energies << control2 << " ";
        energies << energyout << "\n";

        if(i >= stepstart)
        {
            if(i % frequency == 0)
            {
                force_energy_lj(
                                rx,
                                ry,
                                rz,
                                ax,
                                ay,
                                az,
                                l,
                                T,
                                rho,
                                atoms,
                                periodic,
                                dummy,
                                P
                                );

                pressuresum += P;
                count++;

                printf("\n %Lf \n", P);
            }
        }

    }
    P = (long double) pressuresum/count;
    printf("%Lf", P);
    // P *= epsilon/pow(sigma, 3.0);
    pressures << P << "\n";
    energies.close();
    pressures.close();
}
