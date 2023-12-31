#include <omp.h>

void omp_one_regular_time_step (int nx, int ny, double dx, double dy, double dt,
				double **u_new, double **u, double **u_prev)
{
  int i, j, j_left, j_right, i_below, i_above;
  double dx_const, dy_const;

  dx_const = dt*dt/(8*dx*dx);
  dy_const = dt*dt/(8*dy*dy);
  #pragma omp parallel for default(shared) private(j, j_left, j_right, i_below, i_above)\
   firstprivate(dx_const, dy_const)
  for (i = 0; i < ny; i++){
      if (i == 0){
      i_below = 1;
      }
    else{
      i_below = i-1;
    }
    if (i == ny-1){
      i_above = ny-2;     
      }
    else{
      i_above = i + 1;
    }
    for (j = 0; j < nx; j++){
      if (j == 0){
        j_left = 1;
      }
      else{
        j_left = j - 1;
      }
      if (j == nx-1){
        j_right = nx - 2;
      }
      else{
        j_right = j + 1;
      }
      u_new[i][j] = 
      2*u[i][j] + dx_const*(u[i][j_left] - 2*u[i][j] + u[i][j_right])
       + dy_const*(u[i_below][j] -2*u[i][j] + u[i_above][j])-u_prev[i][j]; 
    }
  }
}
