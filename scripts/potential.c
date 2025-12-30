#include <math.h>
#include "../include/config.h"
#include "../include/potential.h"
#include "../include/field.h"
#include "../include/grid.h"

double coulomb_potential(int i, int j, int k) 
{
    int idx = INDEX(i, j, k);
    return V_coulomb[idx];
}

double external_potential(int i, int j, int k, int time) 
{
    return 0;
}

