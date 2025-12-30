#include "../include/config.h"
#include "../include/hamiltonian.h"
#include "../include/timestep.h"
#include <complex.h>

void apply_stencil(complex double* wavefunction, complex double* new_wavefunction, int Nz_local) 
{
    for (int k = 1; k < Nz_local - 1; k++) 
    {
        for (int j = 1; j < Ny - 1; j++) 
        {
            for (int i = 1; i < Nx - 1; i++) 
            {
                int idx = INDEX(i, j, k);
                double H_psi = hamiltonian_operator(i, j, k);
                new_wavefunction[idx] = wavefunction[idx] - I * dt * H_psi * wavefunction[idx];
            }
        }
    }
}

