#include "../include/config.h"
#include "../include/grid.h"
#include <complex.h>

void initialize_grid(complex double* grid, int Nz_local) 
{
    Nz_local += 2; // Account for ghost cells

    for (int k = 0; k < Nz_local; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int idx = INDEX(i, j, k);
                grid[idx] = 0.0 + 0.0 * I;
            }
        }
    }
}

void initialize_r(double* r, int Nz_local) 
{
    Nz_local += 2; // Account for ghost cells

    for (int k = 0; k < Nz_local; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int idx = INDEX(i, j, k);
                double x = i * dx;
                double y = j * dx;
                double z = k * dx;
                r[idx] = sqrt(x * x + y * y + z * z + eps * eps);
            }
        }
    }
}

void initialize_potential(double* V_coulomb, double* r, int Nz_local) 
{
    Nz_local += 2; // Account for ghost cells

    for (int k = 0; k < Nz_local; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int idx = INDEX(i, j, k);
                V_coulomb[idx] = -1.0 / r[idx];
            }
        }
    }
}
