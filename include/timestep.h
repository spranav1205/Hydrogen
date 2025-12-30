#ifndef HYDROGEN_TIMESTEP_H
#define HYDROGEN_TIMESTEP_H

#include <complex.h>

void apply_stencil(complex double* wavefunction, complex double* new_wavefunction, int Nz_local, int time);
void normalize_wavefunction(complex double* wavefunction, int Nz_local);

#endif /* HYDROGEN_TIMESTEP_H */
