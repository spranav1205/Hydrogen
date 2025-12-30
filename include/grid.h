#ifndef HYDROGEN_GRID_H
#define HYDROGEN_GRID_H

#include <complex.h>

#include "config.h"

#define INDEX(i, j, k) ((i) + Nx * ((j) + Ny * (k)))

void initialize_grid(complex double* grid, int Nz_local);
void initialize_r(double* r, int Nz_local);
void initialize_potential(double* V_coulomb, double* r, int Nz_local);

#endif /* HYDROGEN_GRID_H */
