#ifndef HYDROGEN_FIELD_H
#define HYDROGEN_FIELD_H

#include <complex.h>
#include <stdlib.h>
#include "config.h"

extern double _Complex *psi;
extern double _Complex *psi_new;
extern double *r;
extern double *V_coulomb;

void allocate_fields(int Nz_local);
void free_fields();

#endif /* HYDROGEN_FIELD_H */