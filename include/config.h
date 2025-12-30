#ifndef HYDROGEN_CONFIG_H
#define HYDROGEN_CONFIG_H

#define Nx 128
#define Ny 128
#define Nz 128

#define Nt 100

#define kinetic 1.0 /* Factor for kinetic energy: hbar^2 / (2m) */

#define output_interval 10

#define dx 0.01
#define dt 0.001

extern double eps;

#endif /* HYDROGEN_CONFIG_H */