#include "config.h"
#include "grid.h"

#include "complex.h"

void initialize_grid(complex double* grid, int Nz_local) 
{
    for (int k = 0; k < Nz_local; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                int idx = INDEX(i, j, k);
                grid[idx] = 0.0 + 0.0 * I;
            }
        }
    }
}

