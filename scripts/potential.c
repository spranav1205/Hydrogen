#include <math.h>
#include "config.h"
#include "potential.h"
#include "field.h"

double coulomb_potential(int i, int j, int k) 
{
    int idx = INDEX(i, j, k);
    return -1.0/V_coulomb[idx];
}

double external_potential(int i, int j, int k, int time) 
{
    // Placeholder for an external potential, e.g., harmonic oscillator
    int idx = INDEX(i, j, k);

    double k_spring = 0.1; // Spring constant
    return 0.5 * k_spring * r[idx] * r[idx];
}

