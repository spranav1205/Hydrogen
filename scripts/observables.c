#include "config.h"
#include "hamiltonian.h"
#include "observable.h"
#include "field.h"
#include <stddef.h>
#include <mpi.h>
#include <math.h>
#include <complex.h>

MPI_Datatype create_observable_type() {
    MPI_Datatype MPI_OBSERVABLE;
    
    int blocklengths[3] = {1, 3, 1};
    MPI_Aint displacements[3];
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    
    // Calculate displacements
    displacements[0] = offsetof(Observable, probability_density);
    displacements[1] = offsetof(Observable, position);
    displacements[2] = offsetof(Observable, energy);
    
    MPI_Type_create_struct(3, blocklengths, displacements, types, &MPI_OBSERVABLE);
    MPI_Type_commit(&MPI_OBSERVABLE);
    
    return MPI_OBSERVABLE;
}


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

    for (int k = 1; k < Nz_local - 1; k++) 
    {
        for (int j = 1; j < Ny - 1; j++) 
        {
            for (int i = 1; i < Nx - 1; i++) 
            {
                Observable local_obs = measure(wavefunction, i, j, k);
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

Observable gather_observables(int Nz_local, MPI_Datatype MPI_OBSERVABLE)
{
    Observable local_obs, global_obs;
    
    local_obs = expectation(psi, Nz_local);
    MPI_Allreduce(&local_obs, &global_obs, 1, MPI_OBSERVABLE, MPI_SUM, MPI_COMM_WORLD);

    return global_obs;
}

//TODO: Add normalization after gathering observables