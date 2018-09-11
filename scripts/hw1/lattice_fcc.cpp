/*----------------------------------------------------
This script was written by Lane Schultz for MSE 760
Cite me if you got help here
Homework 1
----------------------------------------------------*/

#include <iostream>
using namespace std;

// Return the x, y, and z coordinates of atoms in FCC
float** lattice_fcc(int n, float a)
{
    int initialatoms = 4;                  // The number of attoms in one unit cell
    int dimensions = 3;                    // The dimensions in space
    const int size = 2*n;                  // The size of the arrays
    const int atoms = initialatoms*n*n*n;  // The number of atoms in FCC

    float** matrix = new float*[atoms];    // Matrix to hold the output
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

    // Output
    for(int i = 0; i < atoms; i++)
    {
        matrix[i][0] = x[i];
        matrix[i][1] = y[i];
        matrix[i][2] = z[i];
    }

    return matrix;

    // Delete matrix
    for(size_t i = 0; i < atoms; i++)
        delete matrix[i];
    delete matrix;
}
