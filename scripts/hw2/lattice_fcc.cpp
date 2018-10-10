/*----------------------------------------------------
This script creates an fcc lattice in 3D
----------------------------------------------------*/

// Return the x, y, and z coordinates of atoms in FCC
long double **lattice_fcc(int n, long double a, int *number, int *dim)
{
    int initialatoms = 4;                  // The number of attoms in one unit cell
    int dimensions = 3;                    // The dimensions in space
    const int size = 2*n;                  // The size of the arrays
    const int atoms = initialatoms*n*n*n;  // The number of atoms in FCC

    *number = atoms;
    *dim = dimensions;

    long double **matrix = new long double *[atoms]; // Matrix to hold the output
    for(size_t i = 0; i < atoms; i++)
        matrix[i] = new long double [dimensions];

    // First four atoms
    long double xi[] = {a*0.0, a*0.5, a*0.5, a*0.0};
    long double yi[] = {a*0.0, a*0.5, a*0.0, a*0.5};
    long double zi[] = {a*0.0, a*0.0, a*0.5, a*0.5};

    // The coordinates of the problem
    long double x[atoms] = {};
    long double y[atoms] = {};
    long double z[atoms] = {};

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

    // Delete pointers
    delete number;
    delete dim;

    // Delete matrix
    for(size_t i = 0; i < atoms; i++)
        delete matrix[i];
    delete matrix;
}
