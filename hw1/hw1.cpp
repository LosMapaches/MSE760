/*----------------------------------------------------
This script was written by Lane Schultz for MSE 760
Cite me if you got help here
Homework 1
----------------------------------------------------*/

#include <iostream>
#include <fstream>
using namespace std;

// Function declaration
int lattice_fcc(int n, float a);

// Return the x, y, and z coordinates of atoms in FCC
int lattice_fcc(int n, float a)
{
    std::ofstream outfile;               // Define the output file
    outfile.open("coordinates.txt");     // The output file name
    // Give the output headers
    outfile << "atom x y z\n";

    const int dimensions = 3;            // The dimensions in space
    const int size = 2*n;                // The size of the arrays
    const int atoms = 4*n*n*n;             // The number of atoms in FCC

    float* positions = new float[size];  // Positions along a line
    float** matrix = new float*[atoms];  // Start a matrix to hold values

    // Create N dimensions for each element
    for(size_t i = 0; i < atoms; i++)
	matrix[i] = new float[dimensions];

    // Gather distances in a line
    for(int i = 0; i < size; i++)
	positions[i] = a*i/2.0;

    // First four atoms
    float atom1[] = {a*0, a*0, a*0};
    float atom2[] = {a*0.5, a*0.5, a*0};
    float atom3[] = {a*0, a*0.5, a*0.5};
    float atom4[] = {a*0.5, a*0, a*0.5};

    // Alternate zeroes in coordinates
    for(int i = 0; i < n; i++)
    {

        outfile << i;
        outfile << " ";

	for(int j = 0; j < dimensions; j++)
	{
	    matrix[i][j] = atom1[j]+i*a;
	    outfile << matrix[i][j];
	    outfile << " ";
	}

	outfile << "\n";
        outfile << i;
        outfile << " ";

        for(int j = 0; j < dimensions; j++)
        {
            matrix[i+1][j] = atom2[j]+i*a;
            outfile << matrix[i+1][j];
            outfile << " ";
        }

        outfile << "\n";
        outfile << i;
        outfile << " ";

        for(int j = 0; j < dimensions; j++)
        {
            matrix[i+2][j] = atom3[j]+i*a;
            outfile << matrix[i+2][j];
            outfile << " ";
        }

        outfile << "\n";
        outfile << i;
        outfile << " ";

        for(int j = 0; j < dimensions; j++)
        {
            matrix[i+3][j] = atom4[j]+i*a;
            outfile << matrix[i+3][j];
            outfile << " ";
        }

        outfile << "\n";

    }

    // Close the output file
    outfile.close();

    // Delete matrix to free memory allocation
    delete positions;
    for(size_t i=0; i < size; i++)
        delete matrix[i];
}

int main(void)
{
    int n = 5;                         // Number of unit cells
    float a = 1;                   // Lattice constant 5.256 //4.04092655671750
    float result;

    result = lattice_fcc(n, a); 
}

