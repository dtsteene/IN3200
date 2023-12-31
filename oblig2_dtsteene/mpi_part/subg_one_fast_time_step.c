#include <mpi.h>
void subg_one_fast_time_step(int my_rank, int P, int nx, int my_ny, double dx, double dy, double dt, double **my_u_new,
                             double **my_u, double **my_u_prev)
{
    double dx_const, dy_const;
    int i, j;

    dx_const = dt * dt / (8 * dx * dx);
    dy_const = dt * dt / (8 * dy * dy);

    // No matter the rank we need to calculate the middle
    for (i = 1; i < my_ny; i++)
    {
        // left side
        my_u_new[i][0] = 2 * my_u[i][0] + dx_const * (2 * my_u[i][1] - 2 * my_u[i][0]) +
                         dy_const * (my_u[i - 1][0] - 2 * my_u[i][0] + my_u[i + 1][0]) - my_u_prev[i][0];
        // right side
        my_u_new[i][nx - 1] = 2 * my_u[i][nx - 1] + dx_const * (2 * my_u[i][nx - 2] - 2 * my_u[i][nx - 1]) +
                              dy_const * (my_u[i - 1][nx - 1] - 2 * my_u[i][nx - 1] + my_u[i + 1][nx - 1]) -
                              my_u_prev[i][nx - 1];
        // The middle
        for (j = 1; j < nx - 1; j++)
        {
            my_u_new[i][j] = 2 * my_u[i][j] + dx_const * (my_u[i][j - 1] - 2 * my_u[i][j] + my_u[i][j + 1]) +
                             dy_const * (my_u[i - 1][j] - 2 * my_u[i][j] + my_u[i + 1][j]) - my_u_prev[i][j];
        }
    }

    if (my_rank == 0)
    {
        // Here we also need to calculate the bottom

        // bottom left corner
        my_u_new[0][0] = 2 * my_u[0][0] + dx_const * (2 * my_u[0][1] - 2 * my_u[0][0]) +
                         dy_const * (2 * my_u[1][0] - 2 * my_u[0][0]) - my_u_prev[0][0];
        // bottom right corner
        my_u_new[0][nx - 1] = 2 * my_u[0][nx - 1] + dx_const * (2 * my_u[0][nx - 2] - 2 * my_u[0][nx - 1]) +
                              dy_const * (2 * my_u[1][nx - 1] - 2 * my_u[0][nx - 1]) - my_u_prev[0][nx - 1];
        // bottom side
        for (j = 1; j < nx - 1; j++)
        {
            my_u_new[0][j] = 2 * my_u[0][j] + dx_const * (my_u[0][j - 1] - 2 * my_u[0][j] + my_u[0][j + 1]) +
                             dy_const * (2 * my_u[1][j] - 2 * my_u[0][j]) - my_u_prev[0][j];
        }
    }
    else if (my_rank == P - 1)
    {
        // we are in the top rank and need to calculate the top

        // top left corner
        my_u_new[my_ny][0] = 2 * my_u[my_ny][0] + dx_const * (2 * my_u[my_ny][1] - 2 * my_u[my_ny][0]) +
                             dy_const * (2 * my_u[my_ny - 1][0] - 2 * my_u[my_ny][0]) - my_u_prev[my_ny][0];
        // top right corner
        my_u_new[my_ny][nx - 1] =
            2 * my_u[my_ny][nx - 1] + dx_const * (2 * my_u[my_ny][nx - 2] - 2 * my_u[my_ny][nx - 1]) +
            dy_const * (2 * my_u[my_ny - 1][nx - 1] - 2 * my_u[my_ny][nx - 1]) - my_u_prev[my_ny][nx - 1];

        // top side
        for (j = 1; j < nx - 1; j++)
        {
            my_u_new[my_ny][j] = 2 * my_u[my_ny][j] +
                                 dx_const * (my_u[my_ny][j - 1] - 2 * my_u[my_ny][j] + my_u[my_ny][j + 1]) +
                                 dy_const * (2 * my_u[my_ny - 1][j] - 2 * my_u[my_ny][j]) - my_u_prev[my_ny][j];
        }
    }
    else
    {
        i = my_ny;
        my_u_new[i][0] = 2 * my_u[i][0] + dx_const * (2 * my_u[i][1] - 2 * my_u[i][0]) +
                         dy_const * (my_u[i - 1][0] - 2 * my_u[i][0] + my_u[i + 1][0]) - my_u_prev[i][0];
        // right side
        my_u_new[i][nx - 1] = 2 * my_u[i][nx - 1] + dx_const * (2 * my_u[i][nx - 2] - 2 * my_u[i][nx - 1]) +
                              dy_const * (my_u[i - 1][nx - 1] - 2 * my_u[i][nx - 1] + my_u[i + 1][nx - 1]) -
                              my_u_prev[i][nx - 1];
        // The middle
        for (j = 1; j < nx - 1; j++)
        {
            my_u_new[i][j] = 2 * my_u[i][j] + dx_const * (my_u[i][j - 1] - 2 * my_u[i][j] + my_u[i][j + 1]) +
                             dy_const * (my_u[i - 1][j] - 2 * my_u[i][j] + my_u[i + 1][j]) - my_u_prev[i][j];
        }
    }
}
