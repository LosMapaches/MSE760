/*----------------------------------------------------
This script was written by Lane Schultz for MSE 760
Homework 3
----------------------------------------------------*/

#include <fstream>                         // Use output file
#include <cstdlib>                         // For random numbers
#include <math.h>                          // For powers

#include "lattice_fcc.cpp"                 // FCC coordinates
#include "reduced_units.cpp"               // Reduce units
#include "unreduced_units.cpp"             // Unreduce units
#include "force_energy_lj.cpp"             // Forces and potential energy
#include "mcmove.cpp"                      // Monte Carlo

main()
{
    long double dummy = 0.0;               // Dummy variable to store useless data
    long double m = 6.6e-26;               // Mass [-/atom]
    long double sigma = 3.4e-10;           // Length [m]
    long double epsilon = 0.0104;          // Energy [eV]
    long double k = 8.6173303e-5;          // Boltzmann constant [eV/K]

    int n = 7;                             // Number of units cells
    int atoms = n*n*n*4;                   // Number of atoms
    int steps = 1e6;                       // Number of simulation steps

    // Reduced units
    long double Tred = 2.0;
    long double rhored = 0.40;             // Reduced density
    long double vred = atoms/rhored;
    long double lred = cbrt(vred);
    long double ared = lred/n;

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
    long double P;
    long double pressuresum = 0.0;

    // Atom displacement
    long double delta = 0.2;  // Beginning displacement criterion

    // If move is accepted
    int accept;
    int control1 = 0;  // For the summation of acceptances
    long double control2 = 0.0;  // For determining the percent acceptance

    // Coordinates for FCC lattice
    lattice_fcc(n, ared, rx, ry, rz);

    // Calculate the cohesive energy
    force_energy_lj(
                    rx,
                    ry,
                    rz,
                    ax,
                    ay,
                    az,
                    lred,
                    Tred,
                    rhored,
                    vred,
                    atoms,
                    periodic,
                    cohesive,
                    dummy
                    );

    int frequency = 100;  // The data acquisition rate
    int count = 0;  // The number of times data is taken
    int stepstart = 4e5;  // Approximate number of steps after settling

    // Start Monte Carlo
    srand(time(NULL));  // Generate random seed

    std::ofstream energies;
    energies.open("./0p40rho_energies.txt");

    std::ofstream pressures;
    pressures.open("./0p40rho_pressures.txt");

    energies << "Step AcceptanceRate Energy[eV/atom]\n";
    pressures << "Pressure_virial Pressure_ideal\n";

    for(int i = 1; i <= steps; i++)
    {
       	mcmove(atoms, lred, Tred, delta, rx, ry, rz, cohesive, periodic, energyout, accept);
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
                                lred,
                                Tred,
                                rhored,
                                vred,
                                atoms,
                                periodic,
                                dummy,
                                P
                                );

                pressuresum += P;
                count++;
            }
        }
    }
    P = pressuresum/count;
    printf("Pressure virial: %Lf \n", P);
    printf("Pressure ideal: %Lf \n", rhored*Tred);
    // P *= epsilon/pow(sigma, 3.0);
    pressures << P "\n";
    pressures << rhored*Tred << "\n";
    energies.close();
    pressures.close();
}
