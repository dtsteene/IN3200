void one_regular_time_step (int nx, int ny, double dx, double dy, double dt,
                            double **u_new, double **u, double **u_prev)
{
  int j_left, j_right, i_below, i_above, i, j;
  double dx_const, dy_const;

  dx_const = dt*dt/(8*dx*dx);
  dy_const = dt*dt/(8*dy*dy);
  for (i = 0; i < ny; i++){
    i_below = (i == 0) ? 1 : i-1;
    i_above = (i == ny - 1) ? ny - 2 : i + 1;
    for (j = 0; j < nx; j++){
      j_left = (j == 0) ? 1 : j - 1;
      j_right = (j == nx-1) ? nx - 2: j + 1;
      u_new[i][j] = 
      2*u[i][j] + dx_const*(u[i][j_left] - 2*u[i][j] + u[i][j_right])
       + dy_const*(u[i_below][j] -2*u[i][j] + u[i_above][j])-u_prev[i][j]; 
    }
  }
}
