#ifndef functions_H__
#define functions_H__

void allocate_2D_array (double ***array_ptr, int nx, int ny);

void deallocate_2D_array (double **array);

void first_time_step (int nx, int ny, double dx, double dy, double dt,
		      double **u, double **u_prev);

void one_regular_time_step (int nx, int ny, double dx, double dy, double dt,
			    double **u_new, double **u, double **u_prev);

void one_fast_time_step (int nx, int ny, double dx, double dy, double dt,
			 double **u_new, double **u, double **u_prev);

double compute_numerical_error (int nx, int ny, double dx, double dy,
				double t_value, double **u);

#endif
