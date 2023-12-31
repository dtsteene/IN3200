#define _USE_MATH_DEFINES
#include "../serial_part/functions.h"
#include "parallel_functions.h"
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define Pi M_PI

int main(int nargs, char **args)
{
    int nx = 1001, ny = 1001;
    double T = 2.0; // default values
    int i, j;
    double dx, dy, dt, t;
    double **my_u, **my_u_new, **my_u_prev, **my_tmp_ptr, start_time, finish_time, final_time;
    int my_rank, my_ny, my_offset, P, has_neigh_below, has_neigh_above;

    MPI_Init(&nargs, &args);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    if (my_rank == 0)
    {
        if (nargs > 1) // if a new value of nx is provided on the command line
            nx = ny = atoi(args[1]);
        if (nargs > 2) // if a new value of T is provided on the command line
            T = atof(args[2]);
    }
    // let process 0 broadcast values of nx, ny and T to all other processes
    MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T, 1, MPI_INT, 0, MPI_COMM_WORLD);

    has_neigh_below = (my_rank > 0) ? 1 : 0;
    has_neigh_above = (my_rank < P - 1) ? 1 : 0;
    my_offset = (my_rank * ny) / P;
    my_ny = ((my_rank + 1) * ny) / P - my_offset;
    allocate_2D_array(&my_u, nx, (my_ny + has_neigh_below + has_neigh_above));
    allocate_2D_array(&my_u_new, nx, (my_ny + has_neigh_below + has_neigh_above));
    allocate_2D_array(&my_u_prev, nx, (my_ny + has_neigh_below + has_neigh_above));
    dx = 1.0 / (nx - 1);
    dy = 1.0 / (ny - 1);
    dt = 2.0 * dx; // maximum value allowed for the time step size

    // start timing
    start_time = MPI_Wtime();

    // set initial condition
    for (i = 0; i < my_ny; i++)
        for (j = 0; j < nx; j++)
            my_u_prev[i + has_neigh_below][j] = cos(2.0 * Pi * j * dx) * cos(2.0 * Pi * (i + my_offset) * dy);

    communicate_above_below(my_rank, P, nx, my_ny, my_u_prev);

    subg_first_time_step(my_rank, P, nx, my_ny, dx, dy, dt, my_u, my_u_prev);

    // compute the remaining time steps
    t = dt;
    while (t < T)
    {
        t += dt;
        communicate_above_below(my_rank, P, nx, my_ny, my_u);
        subg_one_fast_time_step(my_rank, P, nx, my_ny, dx, dy, dt, my_u_new, my_u, my_u_prev);
        /* pointer swaps */
        my_tmp_ptr = my_u_prev;
        my_u_prev = my_u;
        my_u = my_u_new;
        my_u_new = my_tmp_ptr;
    }
    printf("my_rank=%d, nx=%d, my_ny=%d, T=%g, dt=%g, error=%e\n", my_rank, nx, my_ny, t, dt,
           all_compute_numerical_error(my_rank, my_offset, P, nx, my_ny, dx, dy, t, my_u));
    // stop timing
    finish_time = MPI_Wtime();
    final_time = finish_time - start_time;
    if (my_rank == 0)
    {
        printf("time %f\n", final_time);
    }
    // deallocate arrays my_u_new, my_u, my_u_prev
    deallocate_2D_array(my_u_new);
    deallocate_2D_array(my_u);
    deallocate_2D_array(my_u_prev);
    MPI_Finalize();
    return 0;
}
