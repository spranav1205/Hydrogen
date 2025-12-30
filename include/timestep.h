#ifndef HYDROGEN_TIMESTEP_H
#define HYDROGEN_TIMESTEP_H

#include <complex.h>

void apply_stencil(complex double* wavefunction, complex double* new_wavefunction, int Nz_local, int time);

#endif /* HYDROGEN_TIMESTEP_H */
