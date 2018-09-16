/*----------------------------------------------------
This script was written by Lane Schultz for MSE 760
Cite me if you got help here
Homework 1
----------------------------------------------------*/

#include <fstream>                         // Use output file
#include <iostream>                        // User input
#include "lattice_fcc.cpp"                 // FCC creator
#include "energy_lj.cpp"                   // LJ energy

main()
{
    std::ofstream outfile;                 // Define the output file
    outfile.open("coordinates.txt");       // The output file name
    std::ofstream outfilecoh;
    outfilecoh.open("cohesive.txt");       // Cohesive energy vs size

    long double a = 5.256e-10;             // Lattice constant [m]
    long double m = 6.6e-23;               // Mass [g]
    long double sigma = 3.4e-10;           // Length [m]
    long double epsilon = 0.0104;          // Energy [eV]


    // Set a loop to spit out cohesive energy for multiple system sizes
    for(int k = 1; k <= 5; k++)
    {
    int n = k;                             // Number of units cells [-]

    /* Grab the atom coordinates for FCC structure */
    int dimensions;                        // Dimension of problem
    int atoms;                             // Number of atoms

    long double **array = lattice_fcc(n, a, &atoms, &dimensions);

    for(int i = 0; i < atoms; i++)
    {
        for(int j = 0; j < dimensions; j++)
        {
            outfile << array[i][j] << " ";
        }
        outfile << "\n";
    }

    outfile.close();                       // Close output file

    long double Utot = energy_lj(array, atoms, sigma, epsilon, m);
    long double Ucoh = Utot/atoms;

    outfilecoh << k << " ";
    outfilecoh << Ucoh << "\n";
    }
    outfilecoh.close();
}

