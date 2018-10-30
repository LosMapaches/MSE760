/*----------------------------------------------------
This script creates an fcc lattice in 3D
----------------------------------------------------*/

// Return the x, y, and z coordinates of atoms in FCC
void lattice_fcc(
                 int n,
                 long double a,
                 long double rx[],
                 long double ry[],
                 long double rz[]
                 )
{
    int initialatoms = 4;                  // The number of attoms in one unit cell
    int dimensions = 3;                    // The dimensions in space
    const int size = 2*n;                  // The size of the arrays
    int atoms = initialatoms*n*n*n;        // The number of atoms

    // First four atoms
    long double rxi[] = {a*0.0, a*0.5, a*0.5, a*0.0};
    long double ryi[] = {a*0.0, a*0.5, a*0.0, a*0.5};
    long double rzi[] = {a*0.0, a*0.0, a*0.5, a*0.5};

    // Clear the coordinates of atoms
    for(int i = 0; i < atoms; i++)
    {
        rx[i] = 0.0;
        ry[i] = 0.0;
        rz[i] = 0.0;
    }

    // Assign the first four atoms to coordinates
    for(int i = 0; i < initialatoms; i++)
    {
        rx[i] = rxi[i];
        ry[i] = ryi[i];
        rz[i] = rzi[i];
    }

    // Repeat along x-axis
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < initialatoms; j++)
        {
            int index = (j+i*initialatoms);
            rx[index] = rx[j]+a*i;
            ry[index] = ry[j];
            rz[index] = rz[j];
        }
    }

    // Repeat along y-axis
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n*initialatoms; j++)
        {
            int index = (j+i*n*initialatoms);
            rx[index] = rx[j];
            ry[index] = ry[j]+a*i;
            rz[index] = rz[j];
        }
    }

    // Repeat along z-axis
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n*n*initialatoms; j++)
        {
            int index = (j+i*n*n*initialatoms);
            rx[index] = rx[j];
            ry[index] = ry[j];
            rz[index] = rz[j]+a*i;
        }
    }

}
