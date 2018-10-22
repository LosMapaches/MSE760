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
    // Energies
    long double pe;

    // Coordinates for FCC lattice
    lattice_fcc(n, a, rx, ry, rz);

    // Velocity coordinates
    velocities(vx, vy, vz, atoms, T);
    kinetic[0] = temperature(atoms, vx, vy, vz);

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
                    pe
                    );

    cohesive[0] = pe;

    // Half values
    long double halfdt = dt/2.0;
    long double halfvx[atoms];
    long double halfvy[atoms];
    long double halfvz[atoms];

    for(int step = 1; step <= steps; step++)
    {
        printf("Step: %i\n", step);

        for(int i = 0; i < atoms; i++)
        {
            halfvx[i] = vx[i]+ax[i]*halfdt;
            halfvy[i] = vy[i]+ay[i]*halfdt;
            halfvz[i] = vz[i]+az[i]*halfdt;

            rx[i] += halfvx[i]*dt;
            ry[i] += halfvy[i]*dt;
            rz[i] += halfvz[i]*dt;
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
                        pe
                        );

        cohesive[step] = pe;

        for(int i = 0; i < atoms; i++)
        {
            vx[i] = halfvx[i]+halfdt*ax[i];
            vy[i] = halfvy[i]+halfdt*ay[i];
            vz[i] = halfvz[i]+halfdt*az[i];
        }

        // The kinetic energy
        kinetic[step] = temperature(atoms, vx, vy, vz);
    }
}
