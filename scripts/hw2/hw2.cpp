/*----------------------------------------------------
This script was written by Lane Schultz for MSE 760
Homework 2
----------------------------------------------------*/

#include <fstream>                         // Use output file
#include "lattice_fcc.cpp"                 // FCC creator
#include "energy_lj.cpp"                   // LJ energy

main()
{
    std::ofstream outfile;                 // Define the output file
    outfile.open("coordinates.txt");       // The output file name

    std::ofstream outfilecoh;
    outfilecoh.open("cohesive.txt");       // Cohesive energy vs size

    long double a = 5.7e-10;               // Lattice constant [m]
    long double m = 6.6e-23;               // Mass [g]

    int n = 5;                             // Number of units cells [-]

    /*Reduce the lattice constant*/
    long double ared = reduced_units(m, 1, a);
    long double l = n*ared;                // Side length [m]

    /* Grab the atom coordinates for FCC structure */
    int dimensions;                        // Dimension of problem
    int atoms;                             // Number of atoms

    long double **acc;                     // The acceleration of atoms
    
    outfile.close();                       // Close output file

    // Set a loop to spit out cohesive energy for multiple system sizes
    outfilecoh << "N E_non-periodic [eV] E_periodic [eV] \n";
    for(int k = 1; k <= 1; k++)
    {
        int n = k;                         // Number of units cells [-]
        long double l = n*ared;            // Side length [m]

        /* Grab the atom coordinates for FCC structure */
        int dimensions;                    // Dimension of problem
        int atoms;                         // Number of atoms

        long double **array = lattice_fcc(n, ared, &atoms, &dimensions);

        long double utotper = energy_lj(array, l, 1, atoms);
        long double ucohper = unreduced_units(m, 2, utotper)/atoms;

        outfilecoh << k << " ";
        outfilecoh << ucohper << "\n";
    }
    outfilecoh.close();
}
