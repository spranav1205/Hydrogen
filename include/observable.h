#ifndef HYDROGEN_OBSERVABLE_H
#define HYDROGEN_OBSERVABLE_H

#include <complex.h>

#include "config.h"
#include "grid.h"

typedef struct {
    double position[3];
    double energy;
} Observable;

Observable measure(complex double* wavefunction, int i, int j, int k);

#endif /* HYDROGEN_OBSERVABLE_H */