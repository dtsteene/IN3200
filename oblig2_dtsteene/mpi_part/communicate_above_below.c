#include <mpi.h>
#include <stdio.h>
void communicate_above_below(int my_rank, int P, int nx, int my_ny, double **my_u)
{
    MPI_Request request1, request2;
    MPI_Status status;
    if ((0 < my_rank) && (my_rank < P - 1))
    {
        MPI_Isend(&my_u[1][0], nx, MPI_DOUBLE, my_rank - 1, 2, MPI_COMM_WORLD, &request1);
        MPI_Irecv(&my_u[my_ny + 1][0], nx, MPI_DOUBLE, my_rank + 1, 2, MPI_COMM_WORLD, &request1);
        MPI_Isend(&my_u[my_ny][0], nx, MPI_DOUBLE, my_rank + 1, 1, MPI_COMM_WORLD, &request2);
        MPI_Irecv(&my_u[0][0], nx, MPI_DOUBLE, my_rank - 1, 1, MPI_COMM_WORLD, &request2);
    }
    else if (my_rank == 0)
    {
        MPI_Irecv(&my_u[my_ny][0], nx, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD, &request1);
        MPI_Isend(&my_u[my_ny - 1][0], nx, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &request2);
    }
    else if (my_rank == P - 1)
    {
        MPI_Isend(&my_u[1][0], nx, MPI_DOUBLE, P - 2, 2, MPI_COMM_WORLD, &request1);
        MPI_Irecv(&my_u[0][0], nx, MPI_DOUBLE, P - 2, 1, MPI_COMM_WORLD, &request2);
    }
    MPI_Wait(&request1, &status);
    MPI_Wait(&request2, &status);
}