#include "../include/config.h"
#include <complex.h>

double complex *psi;
double complex *psi_new;
double *r;
double *V_coulomb;

void allocate_fields(int Nz_local)
{
    size_t N = Nx * Ny * (Nz_local + 2);

    psi = malloc(N * sizeof(double complex));
    psi_new = malloc(N * sizeof(double complex));
    r = malloc(N * sizeof(double));
    V_coulomb = malloc(N * sizeof(double));
}

void free_fields()
{
    free(psi);
    free(psi_new);
    free(r);
    free(V_coulomb);
}