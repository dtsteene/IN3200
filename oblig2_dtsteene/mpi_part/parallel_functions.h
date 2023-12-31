#ifndef parallel_functions_H__
#define parallel_functions_H__

void communicate_above_below (int my_rank, int P, int nx, int my_ny, double **my_u);

void subg_first_time_step (int my_rank, int P, int nx, int my_ny, double dx, double dy, double dt,
double **my_u, double **my_u_prev);

void subg_one_fast_time_step (int my_rank, int P, int nx, int my_ny, double dx, double dy, double dt,
double **my_u_new, double **my_u, double **my_u_prev);

double all_compute_numerical_error (int my_rank, int my_offset, int P, int nx, int my_ny,
double dx, double dy, double t_value, double **my_u);

#endif