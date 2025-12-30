#include "../include/config.h"
#include "../include/potential.h"
#include "../include/hamiltonian.h"
#include "../include/field.h"
#include "../include/grid.h"

double hamiltonian_operator(int i, int j, int k, int time) 
{
    int idx = INDEX(i, j, k);
    double V = coulomb_potential(i, j, k) + external_potential(i, j, k, time);
    double K = -0.5 * kinetic * (
        (psi[INDEX(i+1, j, k)] - 2 * psi[idx] + psi[INDEX(i-1, j, k)]) / (dx * dx) +
        (psi[INDEX(i, j+1, k)] - 2 * psi[idx] + psi[INDEX(i, j-1, k)]) / (dx * dx) +
        (psi[INDEX(i, j, k+1)] - 2 * psi[idx] + psi[INDEX(i, j, k-1)]) / (dx * dx)
    );

    return K + V;
}
    
