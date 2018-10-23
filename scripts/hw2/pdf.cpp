/*----------------------------------------------------
This scripts calculates the pair distribution function
----------------------------------------------------*/

#include <math.h>                          // Do powers

void pdf(
         long double rx[],
         long double ry[],
         long double rz[],
         long double l,
         long double rho,
         int         atoms,
         int         periodic,
         int         bins,
         long double gr[],
         long double dist[]
         )
{
    // Half lengths
    long double poshalf = l/2.0;
    long double neghalf = -poshalf;

    // The difference between vectors
    long double drx;
    long double dry;
    long double drz;

    // The distance between atoms
    long double distance;

    // RDF parameters
    int ig;
    long double delg = l/bins;
    long double norm;

    // Clear RDF values
    for(int i = 0; i < bins; i++)
    {
        gr[i] = 0.0;
        dist[i] = 0.0;
    }

    // Setup the distances
    for(int i = 0; i < bins; i++)
        dist[i] = i*delg;

    for(int i = 0; i < atoms - 1; i++)
    {
        for(int j = i + 1; j < atoms; j++)
        {
            switch(periodic)
            {
            case 0:                     // No periodic boundary

                drx = rx[i]-rx[j];
                dry = ry[i]-ry[j];
                drz = rz[i]-rz[j];

            break;

            case 1:                     // Periodic boundary

                if(rx[i]-rx[j] < neghalf)
                    drx = rx[i]-rx[j]+l;
                else if(rx[i]-rx[j] > poshalf)
                    drx = rx[i]-rx[j]-l;
                else
                    drx = rx[i]-rx[j];

                if(ry[i]-ry[j] < neghalf)
                    dry = ry[i]-ry[j]+l;
                else if(ry[i]-ry[j] > poshalf)
                    dry = ry[i]-ry[j]-l;
                else
                    dry = ry[i]-ry[j];

                if(rz[i]-rz[j] < neghalf)
                    drz = rz[i]-rz[j]+l;
                else if(rz[i]-rz[j] > poshalf)
                    drz = rz[i]-rz[j]-l;
                else
                    drz = rz[i]-rz[j];

            break;
            }
            distance = sqrt(pow(drx, 2.0)+pow(dry, 2.0)+pow(drz, 2.0));

            ig = (int) ceil(distance/delg);
            gr[ig] += 1;
        }
    }
    for(int i = 0; i < bins; i++)
    {
        norm = pow(dist[i]+delg, 3.0)-pow(dist[i], 3.0);
        norm *= (4.0/3.0)*3.14159*rho;
        gr[i] /= norm;
    }
}
