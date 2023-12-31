#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>
#define Pi M_PI

double omp_compute_numerical_error (int nx, int ny, double dx, double dy,
				    double t_value, double **u)
{
  double u_true;
  double sum = 0;
  int i, j;
  #pragma omp parallel for default(shared) private(j, u_true) \
  firstprivate(t_value, dy, dx) reduction(+:sum)
  for (i = 0; i < ny-1; i++){
    for (j = 0; j < nx-1; j++){
      u_true = cos(2*Pi*dx*j)*cos(2*Pi*dy*i)*cos(Pi*t_value);
      sum += pow(u_true-u[i][j],2);
    }
  }
  return sqrt(dx*dy*sum);
  
}

