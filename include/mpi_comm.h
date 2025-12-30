#ifndef MPI_COMM_H
#define MPI_COMM_H

#include <complex.h>

extern int my_rank;
extern int comm_size;

void exchange_ghost_cells(complex double* wavefunction, int Nz_local);

#endif // MPI_COMM_H