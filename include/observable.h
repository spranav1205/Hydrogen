#ifndef HYDROGEN_OBSERVABLE_H
#define HYDROGEN_OBSERVABLE_H

#include <complex.h>
#include <mpi.h>

#include "config.h"
#include "grid.h"

typedef struct {
    double probability_density;
    double position[3];
    double energy;
    double r_squared;
} Observable;

MPI_Datatype create_observable_type(void);
void observable_reduce(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype);
Observable measure(complex double* wavefunction, int i, int j, int k, int time);
Observable expectation(complex double* wavefunction, int Nz_local, int time);
Observable gather_observables(int Nz_local, int time, MPI_Datatype MPI_OBSERVABLE, MPI_Op MPI_OBS_SUM);

#endif /* HYDROGEN_OBSERVABLE_H */