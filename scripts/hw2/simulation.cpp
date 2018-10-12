/*----------------------------------------------------------------------------
Simulate position, velocity, and acceleration updates for a time.
----------------------------------------------------------------------------*/

void simulate(
              int atoms,
              long double rx[],
              long double ry[],
              long double rz[],
              long double vx[],
              long double vy[],
              long double vz[],
              long double ax[],
              long double ay[],
              long double az[],
              long double l,
              long double dt,
              int steps
              )
{
    // Half timestep
    long double halfdt = dt/2.0;

    // The energies for each step calulated
    long double cohesive[steps];
    long double kinetic[steps];
    long double potential[steps];
    long double total[steps];

    // Temperatures
    long double temp[steps];

    for(int step = 1; step <= steps; step++)
    {
        printf("Step: %i\n", step);

        for(int i = 0; i < atoms; i++)
        {
            vx[i] += ax[i]*halfdt;
            vy[i] += ay[i]*halfdt;
            vz[i] += az[i]*halfdt;

            rx[i] += vx[i]*dt;
            ry[i] += vy[i]*dt;
            rz[i] += vz[i]*dt;
        }
        // The acceleration coordinates for each atom
        long double energy_cohesive = force_energy_lj(
                                                      rx,
                                                      ry,
                                                      rz,
                                                      ax,
                                                      ay,
                                                      az,
                                                      l,
                                                      atoms,
                                                      1
                                                      );

    }
}
