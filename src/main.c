#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <mpi.h>
#include "../include/observable.h"
#include "../include/hamiltonian.h"
#include "../include/grid.h"
#include "../include/config.h"
#include "../include/field.h"
#include "../include/mpi_comm.h"
#include "../include/timestep.h"

int my_rank, comm_size;

int main()
{
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    MPI_Op MPI_OBS_SUM;
    MPI_Op_create(observable_reduce, 1, &MPI_OBS_SUM);
    MPI_Datatype MPI_OBSERVABLE = create_observable_type();

    int Nz_local = Nz / comm_size; // Simple division of grid in z-direction

    allocate_fields(Nz_local);
    initialize_r(r, Nz_local);
    initialize_grid(psi, r, Nz_local);
    initialize_potential(V_coulomb, r, Nz_local);
    normalize_wavefunction(psi, Nz_local); 

    // Report initial state before any time stepping
    Observable obs0 = gather_observables(Nz_local, 0, MPI_OBSERVABLE, MPI_OBS_SUM);
    if (my_rank == 0) {
        double energy_eV0 = obs0.energy * 27.2;
        printf("Time step %d: <E> = %.6f a.u. = %.3f eV\n", 0, obs0.energy, energy_eV0);
        printf("               <x> = %.4f, <y> = %.4f, <z> = %.4f Bohr\n", 
               obs0.position[0], obs0.position[1], obs0.position[2]);
        printf("               <r> = %.6f Bohr\n", sqrt(obs0.r_squared));
        printf("               Expected (hydrogen 1s): E₀ = -0.5 a.u. = -13.6 eV, <r> = 1.5 Bohr\n");
    }

    for (int time=0; time<Nt; time++) {
        apply_stencil(psi, psi_new, Nz_local, time);
        exchange_ghost_cells(psi_new, Nz_local);

        // Swap pointers
        complex double* temp = psi;
        psi = psi_new;
        psi_new = temp;

        normalize_wavefunction(psi, Nz_local);

        int step = time + 1; // after one update
        if (step % output_interval == 0) {
            Observable obs = gather_observables(Nz_local, step, MPI_OBSERVABLE, MPI_OBS_SUM);
            if (my_rank == 0) {
                // Convert atomic units to eV for display
                // 1 Hartree = 27.2 eV
                double energy_eV = obs.energy * 27.2;
                printf("Time step %d: <E> = %.6f a.u. = %.3f eV\n", step, obs.energy, energy_eV);
                printf("               <x> = %.4f, <y> = %.4f, <z> = %.4f Bohr\n", 
                       obs.position[0], obs.position[1], obs.position[2]);
                printf("               <r> = %.6f Bohr\n", sqrt(obs.r_squared));
            }
        }
    }
    
    free_fields();
    MPI_Type_free(&MPI_OBSERVABLE);
    MPI_Op_free(&MPI_OBS_SUM);
    MPI_Finalize();
    return 0;
}