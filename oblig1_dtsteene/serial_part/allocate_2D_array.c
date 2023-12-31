#include <stdlib.h>

void allocate_2D_array (double ***array_ptr, int nx, int ny)
{
  /*
  double* temp_storage = (double*)malloc(nx*ny*sizeof(double));
  double** temp = (double**)malloc(ny*sizeof(double*));
  for (int i = 0; i < ny; i++)
    temp[i] = &(temp_storage[i*nx]);

  *array_ptr = temp;
  */

  double **temp_arr;
  temp_arr = (double **)malloc(ny * sizeof(double *));
  temp_arr[0] = (double *)malloc(ny * nx * sizeof(double));

  for (int i = 1; i < ny; i ++) {
      temp_arr[i] = &(temp_arr[0][nx * i]);
  }
  *array_ptr = temp_arr;
}
