#define _USE_MATH_DEFINES
#include <math.h>
#include <mpi.h>
#define Pi M_PI

double all_compute_numerical_error(int my_rank, int my_offset, int P, int nx, int my_ny, double dx, double dy,
                                   double t_value, double **my_u)
{
    double u_true;
    double sum = 0, global_sum = 0;
    int has_neigh_below = (my_rank > 0) ? 1 : 0;

    int i, j;
    for (i = 0; i < my_ny - 1; i++)
    {
        for (j = 0; j < nx - 1; j++)
        {
            u_true = cos(2 * Pi * dx * j) * cos(2 * Pi * dy * (i + my_offset)) * cos(Pi * t_value);
            sum += pow(u_true - my_u[i + has_neigh_below][j], 2);
        }
    }
    MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sqrt(dx * dy * global_sum);
}