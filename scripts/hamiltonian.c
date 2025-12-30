#include "../include/config.h"
#include "../include/potential.h"

double hamiltonian_operator(int i, int j, int k, int time) 
{
    int idx = INDEX(i, j, k);
    double V = coulomb_potential(i, j, k) + external_potential(i, j, k, time);
    double K = -0.5 * kinetic * (
        (INDEX(i+1, j, k) - 2 * idx + INDEX(i-1, j, k)) / (dx * dx) +
        (INDEX(i, j+1, k) - 2 * idx + INDEX(i, j-1, k)) / (dx * dx) +
        (INDEX(i, j, k+1) - 2 * idx + INDEX(i, j, k-1)) / (dx * dx)
    );

    return K + V;
}
    
