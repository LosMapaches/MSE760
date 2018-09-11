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

    int initialatoms = 4;          // The number of attoms in one unit cell
    int dimensions = 3;            // The dimensions in space
    const int size = 2*n;                // The size of the arrays
    const int atoms = initialatoms*n*n*n;             // The number of atoms in FCC

    float** matrix = new float*[atoms];  // Start a matrix to hold values

    // Create N dimensions for each element
    for(size_t i = 0; i < atoms; i++)
	matrix[i] = new float[dimensions];

    // First four atoms
    float xi[] = {a*0, a*0.5, a*0.5, a*0};
    float yi[] = {a*0, a*0.5, a*0, a*0.5};
    float zi[] = {a*0, a*0, a*0.5, a*0.5};

    // The coordinates of the problem
    float x[atoms] = {};
    float y[atoms] = {};
    float z[atoms] = {};

    // Initial structure to be repeated
    for(int i = 0; i < initialatoms; i++)
    {
        x[i] = xi[i];
        y[i] = yi[i];
        z[i] = zi[i];
    }

    // Repeat along x-axis
    for(int i = 0; i < n; i++)
    {
    for(int j = 0; j < initialatoms; j++)
    {
        int index = (j+i*initialatoms);
        x[index] = x[j]+a*i;
        y[index] = y[j];
        z[index] = z[j];
    }
    }

    // Repeat along y-axis
    for(int i = 0; i < n; i++)
    {
    for(int j = 0; j < n*initialatoms; j++)
    {
        int index = (j+i*n*initialatoms);
        x[index] = x[j];
        y[index] = y[j]+a*i;
        z[index] = z[j];
    }
    }

    // Repeat along z-axis
    for(int i = 0; i < n; i++)
    {
    for(int j = 0; j < n*n*initialatoms; j++)
    {
        int index = (j+i*n*n*initialatoms);
        x[index] = x[j];
        y[index] = y[j];
        z[index] = z[j]+a*i;
    }
    }

    // Output text file
    for(int i = 0; i < atoms; i++)
    {
        outfile << x[i] << " ";
        outfile << y[i] << " ";
        outfile << z[i] << " ";
        outfile << "\n";
    }

    // Close the output file
    outfile.close();

    // Delete matrix to free memory allocation
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

