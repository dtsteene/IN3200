#include <stdlib.h>

void deallocate_2D_array (double **array)
{
  free(array[0]);
  free(array);
}
