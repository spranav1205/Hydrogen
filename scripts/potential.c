#include <math.h>
#include "config.h"

double coulomb_potential(int i, int j, int k) 
{
    double x = i * dx;
    double y = j * dx;
    double z = k * dx;

    double r = sqrt(x * x + y * y + z * z + eps * eps);

    return 1.0 / r;
}

