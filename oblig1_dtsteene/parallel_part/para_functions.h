#ifndef para_functions_H__
#define para_functions_H__

void omp_first_time_step (int nx, int ny, double dx, double dy, double dt,
			  double **u, double **u_prev);

void omp_one_regular_time_step (int nx, int ny, double dx, double dy, double dt,
				double **u_new, double **u, double **u_prev);

void omp_one_fast_time_step (int nx, int ny, double dx, double dy, double dt,
			     double **u_new, double **u, double **u_prev);

double omp_compute_numerical_error (int nx, int ny, double dx, double dy,
				    double t_value, double **u);


#endif
