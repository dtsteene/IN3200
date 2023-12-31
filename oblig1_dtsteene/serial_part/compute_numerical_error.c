#define _USE_MATH_DEFINES
#include <math.h>

double compute_numerical_error (int nx, int ny, double dx, double dy,
                                double t_value, double **u)
{
  int i, j;
  double u_true;
  double sum = 0;
  for (i = 0; i < ny-1; i++){
    for (j = 0; j < nx-1; j++){
      u_true = cos(2*M_PI*dx*j)*cos(2*M_PI*dy*i)*cos(M_PI*t_value);
      sum += pow(u_true-u[i][j],2);
    }
  }
  return sqrt(dx*dy*sum);
}

