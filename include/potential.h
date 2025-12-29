#ifndef HYDROGEN_POTENTIAL_H
#define HYDROGEN_POTENTIAL_H

#include "config.h"

double coulomb_potential(int i, int j, int k);
double external_potential(int i, int j, int k, int time);

#endif /* HYDROGEN_POTENTIAL_H */

