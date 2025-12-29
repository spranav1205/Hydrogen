#include "config.h"
#include "hamiltonian.h"
#include "observable.h"
#include <math.h>
#include <complex.h>

Observable measure(complex double* wavefunction, int i, int j, int k)
{
    int idx = INDEX(i, j, k);
    double E = hamiltonian_operator(i, j, k);
    double probability_density = creal(wavefunction[idx] * conj(wavefunction[idx]));

    Observable obs;
    obs.probability_density = probability_density;
    obs.position[0] = i * dx;
    obs.position[1] = j * dx;
    obs.position[2] = k * dx;
    obs.energy = E;

    return obs;
}

Observable expectation(complex double* wavefunction, int Nz_local)
{
    Observable obs;
    obs.probability_density = 0.0;
    obs.position[0] = 0.0;
    obs.position[1] = 0.0;
    obs.position[2] = 0.0;
    obs.energy = 0.0;

    double total_probability = 0.0;

    for (int k = 1; k < Nz_local - 1; k++) 
    {
        for (int j = 1; j < Ny - 1; j++) 
        {
            for (int i = 1; i < Nx - 1; i++) 
            {
                Observable local_obs = measure(wavefunction, i, j, k);
                total_probability += local_obs.probability_density;
                obs.probability_density += local_obs.probability_density;
                obs.position[0] += local_obs.position[0] * local_obs.probability_density;
                obs.position[1] += local_obs.position[1] * local_obs.probability_density;
                obs.position[2] += local_obs.position[2] * local_obs.probability_density;
                obs.energy += local_obs.energy * local_obs.probability_density;
            }
        }
    }
    
    return obs;
}