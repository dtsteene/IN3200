void first_time_step (int nx, int ny, double dx, double dy, double dt,
                      double **u, double **u_prev)
{
  int j_left, j_right, i_below, i_above, i, j;
  for (i = 0; i < ny; i++){
    i_below = (i == 0) ? 1 : i-1;
    i_above = (i == ny - 1) ? ny - 2 : i + 1;
    for (j = 0; j < nx; j++){
      j_left = (j == 0) ? 1 : j - 1;
      j_right = (j == nx-1) ? nx - 2: j + 1;
      u[i][j] = u_prev[i][j] + dt*dt/(16*dx*dx)*(u_prev[i][j_left] - 2*u_prev[i][j] + u_prev[i][j_right]) + dt*dt/(16*dy*dy)*(u_prev[i_below][j] -2*u_prev[i][j] + u_prev[i_above][j]); 
    }
  }
}
