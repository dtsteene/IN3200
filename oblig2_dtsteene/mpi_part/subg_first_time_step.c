#include <mpi.h>

void subg_first_time_step(int my_rank, int P, int nx, int my_ny, double dx, double dy, double dt, double **my_u,
                          double **my_u_prev)
{

    int j_left, j_right, i_below, i_above, i, j;
    double dx_const, dy_const;
    dx_const = dt * dt / (16 * dx * dx);
    dy_const = dt * dt / (16 * dy * dy);

    if (my_rank == 0)
    {
        for (i = 0; i < my_ny; i++)
        {
            i_below = (i == 0) ? 1 : i - 1;
            i_above = i + 1;
            for (j = 0; j < nx; j++)
            {
                j_left = (j == 0) ? 1 : j - 1;
                j_right = (j == nx - 1) ? nx - 2 : j + 1;
                my_u[i][j] = my_u_prev[i][j] +
                             dx_const * (my_u_prev[i][j_left] - 2 * my_u_prev[i][j] + my_u_prev[i][j_right]) +
                             dy_const * (my_u_prev[i_below][j] - 2 * my_u_prev[i][j] + my_u_prev[i_above][j]);
            }
        }
    }
    else if (my_rank == P - 1)
    {
        for (i = 1; i < my_ny + 1; i++)
        {
            i_below = i - 1;
            i_above = (i == my_ny) ? my_ny - 1 : i + 1;
            for (j = 0; j < nx; j++)
            {
                j_left = (j == 0) ? 1 : j - 1;
                j_right = (j == nx - 1) ? nx - 2 : j + 1;
                my_u[i][j] = my_u_prev[i][j] +
                             dx_const * (my_u_prev[i][j_left] - 2 * my_u_prev[i][j] + my_u_prev[i][j_right]) +
                             dy_const * (my_u_prev[i_below][j] - 2 * my_u_prev[i][j] + my_u_prev[i_above][j]);
            }
        }
    }
    else
    {
        for (i = 1; i < my_ny + 1; i++)
        {
            i_below = i - 1;
            i_above = i + 1;
            for (j = 0; j < nx; j++)
            {
                j_left = (j == 0) ? 1 : j - 1;
                j_right = (j == nx - 1) ? nx - 2 : j + 1;
                my_u[i][j] = my_u_prev[i][j] +
                             dx_const * (my_u_prev[i][j_left] - 2 * my_u_prev[i][j] + my_u_prev[i][j_right]) +
                             dy_const * (my_u_prev[i_below][j] - 2 * my_u_prev[i][j] + my_u_prev[i_above][j]);
            }
        }
    }
}
