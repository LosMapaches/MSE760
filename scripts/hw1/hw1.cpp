/*----------------------------------------------------
This script was written by Lane Schultz for MSE 760
Cite me if you got help here
Homework 1
----------------------------------------------------*/

#include <iostream>
#include <fstream>
#include "lattice_fcc.cpp"
using namespace std;

int main(void)
{
    std::ofstream outfile;                 // Define the output file
    outfile.open("coordinates.txt");       // The output file name

    int n = 5;                             // Number of unit cells
    int dimensions = 3;                    // Dimension of problem
    int atoms = n*n*n*4;                   // Number of atoms
    float a = 5.256;                       // Lattice constant 5.256 [A]

    float** result = lattice_fcc(n, a);

    for(int i = 0; i < atoms; i++)
    {
        for(int j = 0; j < dimensions; j++)
        {
            outfile << result[i][j] << " ";
        }
        outfile << "\n";
    }

    outfile.close();
}

