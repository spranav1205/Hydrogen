#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <mpi.h>
#include "../include/observable.h"
#include "../include/hamiltonian.h"
#include "../include/config.h"
#include "../include/field.h"

int main()
{
    int my_rank, comm_size;
    
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    

}