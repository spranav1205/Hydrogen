#ifndef HYDROGEN_CONFIG_H
#define HYDROGEN_CONFIG_H

// ============================================================================
// ATOMIC UNITS: ℏ=1, m_e=1, e=1, 4πε₀=1
// Length in Bohr radii (a₀ = 0.529 Å)
// Energy in Hartree (E_h = 2 × 13.6 eV = 27.2 eV)
// Time in ℏ/E_h = 2.419 × 10^-17 s (atomic time unit)
// ============================================================================

// Grid parameters (refined)
#define Nx 500
#define Ny 500
#define Nz 500

#define x_center (Nx - 1) * dx / 2.0
#define y_center (Ny - 1) * dx / 2.0
#define z_center (Nz - 1) * dx / 2.0

#define Nt 100

// Physical parameters in ATOMIC UNITS
#define dx 0.02                    // 0.1 Bohr radii (~0.0529 Å) finer grid
#define dt 0.00005                // smaller dt for stability on finer grid (explicit Euler)
#define dV (dx * dx * dx)        // Volume element

// In atomic units: kinetic = ℏ²/(2m) = 1/2
#define kinetic 0.5

#define output_interval 10

extern double eps;

#endif /* HYDROGEN_CONFIG_H */