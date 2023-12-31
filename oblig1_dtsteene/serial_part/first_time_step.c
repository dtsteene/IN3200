void first_time_step (int nx, int ny, double dx, double dy, double dt,
                      double **u, double **u_prev)
{
  int j_left, j_right, i_below, i_above, i, j;
  for (i = 0; i < ny; i++){
    if (i == 0){
      i_below = 1;
    }
    else{
      i_below = i - 1;
    }
    if (i == ny-1){
      i_above = ny - 2;
    }
    else{
      i_above = i + 1;
    }
    for (j = 0; j < nx; j++){
      if (j == 0){
        j_left = 1;
      }
      else{
        j_left = j-1;
      }
      if (j == nx-1){
        j_right = nx-2;     
        }
      else{
        j_right = j + 1;
      }
      u[i][j] = u_prev[i][j] + dt*dt/(16*dx*dx)*(u_prev[i][j_left] - 2*u_prev[i][j] + u_prev[i][j_right]) + dt*dt/(16*dy*dy)*(u_prev[i_below][j] -2*u_prev[i][j] + u_prev[i_above][j]); 
    }
  }
}
