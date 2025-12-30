#define _DEFAULT_SOURCE
#include "../include/config.h"
#include "../include/grid.h"
#include <math.h>
#include <complex.h>

void initialize_grid(complex double* grid, double* r, int Nz_local) 
{
    for (int k = 0; k < Nz_local+2; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int idx = INDEX(i, j, k);
                
                // 1s orbital: ψ_{1s}(r) = (1/√π) * exp(-r) (in atomic units)
                double psi_1s = (1.0 / sqrt(M_PI)) * exp(-r[idx]);
                grid[idx] = psi_1s + 0.0 * I;
            }
        }
    }
}

void initialize_r(double* r, int Nz_local) 
{
    int z_offset = my_rank * Nz_local; // account for leading ghost cell

    for (int k = 0; k < Nz_local+2; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int idx = INDEX(i, j, k);
                double x = i * dx - x_center;
                double y = j * dx - y_center;
                // IMPORTANT: global z coordinate
                double z = (z_offset + k - 1) * dx - z_center;
                r[idx] = sqrt(x * x + y * y + z * z + eps * eps);
            }
        }
    }
}

void initialize_potential(double* V_coulomb, double* r, int Nz_local) 
{
    for (int k = 0; k < Nz_local+2; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int idx = INDEX(i, j, k);
                // Coulomb potential in atomic units: V(r) = -1/r
                V_coulomb[idx] = -1.0 / r[idx];
            }
        }
    }
}
