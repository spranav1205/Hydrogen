#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
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

    int Nz_local = Nz / comm_size; // Simple division of grid in z-direction

    allocate_fields(Nz_local);
    initialize_grid(psi, Nz_local);
    initialize_r(r, Nz_local);
    initialize_potential(V_coulomb, r, Nz_local);

    for (int time=0; time<Nt; time++) {
        apply_stencil(psi, psi_new, Nz_local, time);
        exchange_ghost_cells(psi_new, Nz_local);

        // Swap pointers
        complex double* temp = psi;
        psi = psi_new;
        psi_new = temp;

        if (time % output_interval == 0) {
            Observable obs = expectation(psi, Nz_local, time);
            if (my_rank == 0) {
                printf("Time step %d: <E> = %f\n", time, obs.energy);
                printf("               <x> = %f, <y> = %f, <z> = %f\n", obs.position[0], obs.position[1], obs.position[2]);
            }
        }
    }
    
    free_fields();
    MPI_Finalize();
    return 0;
}