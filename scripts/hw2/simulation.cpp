/*----------------------------------------------------------------------------
Simulate position, velocity, and acceleration updates for a time.
----------------------------------------------------------------------------*/
#include <stdio.h>

long double **simulate(
                       int atoms,
                       int steps,
                       long double l,
                       long double Ti,
                       long double timestep,
                       long double **r
                       )
{
    // The energies for each step calulated
    long double energy0;
    long double cohesive[steps];
    long double kinetic[atoms];
    long double potential[atoms];
    long double total[atoms];

    // Initial values
    long double checktemp;
    long double **v = velocities(atoms, Ti, &checktemp);
    long double **a = force_energy_lj(r, l, 1, atoms, &energy0);
    cohesive[0] = energy0;

    for(int step = 1; step <= steps; step++)
    {
        printf("Step: %i\n", step);

        for(int i = 0; i < atoms; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                v[i][j] += a[i][j]*timestep/2.0;
                r[i][j] += v[i][j]*timestep;
            }
        }

        a = force_energy_lj(r, l, 1, atoms, &energy0);
        cohesive[step] = energy0;
    }

    return r;

}
