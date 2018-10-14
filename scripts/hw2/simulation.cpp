/*----------------------------------------------------------------------------
Simulate position, velocity, and acceleration updates for a time.
----------------------------------------------------------------------------*/

void simulate(
              int         atoms,
              int         n,
              long double a,
              long double T,
              long double cohesive[],
              long double kinetic[],
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
              int         steps
              )
{

    // Coordinates for FCC lattice
    lattice_fcc(n, a, rx, ry, rz);

    // Velocity coordinates
    long double tempcheck = velocities(vx, vy, vz, atoms, T);

    // The acceleration coordinates for each atom
    force_energy_lj(
                    rx,
                    ry,
                    rz,
                    ax,
                    ay,
                    az,
                    l,
                    atoms,
                    1,
                    cohesive[0]
                    );

    // Half timestep
    long double halfdt = dt/2.0;

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
        force_energy_lj(
                        rx,
                        ry,
                        rz,
                        ax,
                        ay,
                        az,
                        l,
                        atoms,
                        1,
                        cohesive[step]
                        );
    }
}
