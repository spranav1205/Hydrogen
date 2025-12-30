#include "../include/config.h"
#include "../include/hamiltonian.h"
#include "../include/timestep.h"
#include "../include/grid.h"
#include <complex.h>
#include <math.h>
#include <mpi.h>

void apply_stencil(complex double* wavefunction, complex double* new_wavefunction, int Nz_local, int time) 
{
    for (int k = 1; k <= Nz_local; k++) // Skip ghost cells (IMPORTANT)
    {
        for (int j = 1; j < Ny - 1; j++) 
        {
            for (int i = 1; i < Nx - 1; i++) 
            {
                int idx = INDEX(i, j, k);
                complex double H_psi = hamiltonian_operator(i, j, k, time);
                new_wavefunction[idx] = wavefunction[idx] - I * dt * H_psi;
            }
        }
    }
}

void normalize_wavefunction(complex double* wavefunction, int Nz_local)
{
    double local_norm = 0.0;
    for (int k = 1; k <= Nz_local; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int idx = INDEX(i, j, k);
                local_norm += creal(wavefunction[idx] * conj(wavefunction[idx])) * dV; // probability density * dV
            }
        }
    }

    double global_norm = 0.0;
    MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (global_norm <= 0.0) {
        return;
    }
    double scale = 1.0 / sqrt(global_norm);

    // Scale entire array including ghosts to keep boundaries consistent
    int N = Nx * Ny * (Nz_local + 2);
    for (int idx = 0; idx < N; idx++) {
        wavefunction[idx] *= scale;
    }
}

