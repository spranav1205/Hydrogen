#include "../include/config.h"
#include "../include/mpi_comm.h"
#include "../include/field.h"
#include <mpi.h>

void exchange_ghost_cells(complex double* wavefunction, int Nz_local) 
{
    MPI_Request requests[4];
    MPI_Status statuses[4];
    int num_requests = 0; // Waits for these many requests, depending on rank

    if(my_rank > 0 && my_rank < comm_size - 1) 
    {
        // Non-blocking send to previous rank
        // Non-blocking receive from previous rank
        MPI_Isend(&wavefunction[INDEX(0, 0, 1)], Nx * Ny, MPI_C_DOUBLE_COMPLEX, my_rank - 1, 0, MPI_COMM_WORLD, &requests[num_requests++]);
        MPI_Irecv(&wavefunction[INDEX(0, 0, 0)], Nx * Ny, MPI_C_DOUBLE_COMPLEX, my_rank - 1, 1, MPI_COMM_WORLD, &requests[num_requests++]);

        // Non-blocking send to next rank
        // Non-blocking receive from next rank
        MPI_Isend(&wavefunction[INDEX(0, 0, Nz_local)], Nx * Ny, MPI_C_DOUBLE_COMPLEX, my_rank + 1, 1, MPI_COMM_WORLD, &requests[num_requests++]);
        MPI_Irecv(&wavefunction[INDEX(0, 0, Nz_local+1)], Nx * Ny, MPI_C_DOUBLE_COMPLEX, my_rank + 1, 0, MPI_COMM_WORLD, &requests[num_requests++]);
        
        MPI_Waitall(num_requests, requests, statuses);
    }
    else if(my_rank == 0)
    {

        MPI_Isend(&wavefunction[INDEX(0, 0, Nz_local)], Nx * Ny, MPI_C_DOUBLE_COMPLEX, 1, 1, MPI_COMM_WORLD, &requests[num_requests++]);
        MPI_Irecv(&wavefunction[INDEX(0, 0, Nz_local+1)], Nx * Ny, MPI_C_DOUBLE_COMPLEX, 1, 0, MPI_COMM_WORLD, &requests[num_requests++]);

        MPI_Waitall(num_requests, requests, statuses);
    }
    else if(my_rank == comm_size - 1)
    {
        MPI_Isend(&wavefunction[INDEX(0, 0, 1)], Nx * Ny, MPI_C_DOUBLE_COMPLEX, my_rank - 1, 1, MPI_COMM_WORLD, &requests[num_requests++]);
        MPI_Irecv(&wavefunction[INDEX(0, 0, 0)], Nx * Ny, MPI_C_DOUBLE_COMPLEX, my_rank - 1, 0, MPI_COMM_WORLD, &requests[num_requests++]);

        MPI_Waitall(num_requests, requests, statuses);
    }
}