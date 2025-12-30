#include "../include/config.h"
#include "../include/hamiltonian.h"
#include "../include/observable.h"
#include "../include/field.h"
#include <stddef.h>
#include <mpi.h>
#include <math.h>
#include <complex.h>

extern int my_rank;
extern int comm_size;

// Custom reduction function for Observable struct
void observable_reduce(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
{
    Observable *in = (Observable *)invec;
    Observable *inout = (Observable *)inoutvec;
    
    for (int i = 0; i < *len; i++) {
        inout[i].probability_density += in[i].probability_density;
        inout[i].position[0] += in[i].position[0];
        inout[i].position[1] += in[i].position[1];
        inout[i].position[2] += in[i].position[2];
        inout[i].energy += in[i].energy;
        inout[i].r_squared += in[i].r_squared;
    }
}

MPI_Datatype create_observable_type() {
    MPI_Datatype MPI_OBSERVABLE;
    
    int blocklengths[4] = {1, 3, 1, 1};
    MPI_Aint displacements[4];
    MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    
    // Calculate displacements
    displacements[0] = offsetof(Observable, probability_density);
    displacements[1] = offsetof(Observable, position);
    displacements[2] = offsetof(Observable, energy);
    displacements[3] = offsetof(Observable, r_squared);
    
    MPI_Type_create_struct(4, blocklengths, displacements, types, &MPI_OBSERVABLE);
    MPI_Type_commit(&MPI_OBSERVABLE);
    
    return MPI_OBSERVABLE;
}


Observable measure(complex double* wavefunction, int i, int j, int k, int time)
{
    int idx = INDEX(i, j, k);
    complex double H_psi = hamiltonian_operator(i, j, k, time);
    // Energy expectation: <ψ|H|ψ> / <ψ|ψ>
    double E = creal(conj(wavefunction[idx]) * H_psi);
    double probability_density = creal(wavefunction[idx] * conj(wavefunction[idx]));

    // Centered coordinates (global) so <r> should be ~0 for symmetric states
    double x_center = (Nx - 1) * dx / 2.0;
    double y_center = (Ny - 1) * dx / 2.0;
    int Nz_global = Nz / comm_size; 
    int z_global = my_rank * Nz_global + (k - 1); // k runs 1..Nz_local (physical cells)
    double z_center = (Nz - 1) * dx / 2.0;

    Observable obs;
    obs.probability_density = probability_density;
    obs.position[0] = i * dx - x_center;
    obs.position[1] = j * dx - y_center;
    obs.position[2] = z_global * dx - z_center;
    obs.energy = E;
    obs.r_squared = obs.position[0]*obs.position[0] + obs.position[1]*obs.position[1] + obs.position[2]*obs.position[2];

    return obs;
}

Observable expectation(complex double* wavefunction, int Nz_local, int time)
{
    const double dV = dx * dx * dx;
    Observable obs;
    obs.probability_density = 0.0;
    obs.position[0] = 0.0;
    obs.position[1] = 0.0;
    obs.position[2] = 0.0;
    obs.energy = 0.0;
    obs.r_squared = 0.0;

    for (int k = 1; k <= Nz_local; k++) 
    {
        for (int j = 1; j < Ny - 1; j++) 
        {
            for (int i = 1; i < Nx - 1; i++) 
            {
                Observable local_obs = measure(wavefunction, i, j, k, time);
                double weight = local_obs.probability_density * dV;
                obs.probability_density += weight;
                obs.position[0] += local_obs.position[0] * weight;
                obs.position[1] += local_obs.position[1] * weight;
                obs.position[2] += local_obs.position[2] * weight;
                obs.energy += local_obs.energy * weight;
                obs.r_squared += local_obs.r_squared * weight;
            }
        }
    }

    return obs;
}

Observable gather_observables(int Nz_local, int time, MPI_Datatype MPI_OBSERVABLE)
{
    Observable local_obs, global_obs;
    
    local_obs = expectation(psi, Nz_local, time);
    
    // Create custom MPI operation for Observable reduction
    MPI_Op MPI_OBS_SUM;
    MPI_Op_create(observable_reduce, 1, &MPI_OBS_SUM);
    
    MPI_Allreduce(&local_obs, &global_obs, 1, MPI_OBSERVABLE, MPI_OBS_SUM, MPI_COMM_WORLD);
    
    MPI_Op_free(&MPI_OBS_SUM);

    if (global_obs.probability_density > 0) {
        global_obs.position[0] /= global_obs.probability_density;
        global_obs.position[1] /= global_obs.probability_density;
        global_obs.position[2] /= global_obs.probability_density;
        global_obs.energy /= global_obs.probability_density;
        global_obs.r_squared /= global_obs.probability_density;
    }

    return global_obs;
}
