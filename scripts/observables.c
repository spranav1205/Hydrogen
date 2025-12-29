#include "config.h"
#include "observable.h"
#include <math.h>
#include <complex.h>

Observable measure(complex double* wavefunction, int i, int j, int k)
{
    int idx = INDEX(i, j, k);
    double probability_density = creal(wavefunction[idx] * conj(wavefunction[idx]));
    Observable obs;
    obs.position[0] = i * dx;
    obs.position[1] = j * dx;
    obs.position[2] = k * dx;
}