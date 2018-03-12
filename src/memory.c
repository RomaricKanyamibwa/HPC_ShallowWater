#include <stdlib.h>
#include <shalw.h>

void alloc(void) {
  if(my_rank==0)
  {
    hFil = (double *) calloc(2*size_x*size_y, sizeof(double));
    uFil = (double *) calloc(2*size_x*size_y, sizeof(double));
    vFil = (double *) calloc(2*size_x*size_y, sizeof(double));
    hPhy = (double *) calloc(2*size_x*size_y, sizeof(double));
    uPhy = (double *) calloc(2*size_x*size_y, sizeof(double));
    vPhy = (double *) calloc(2*size_x*size_y, sizeof(double));
  }
    hFil_local = (double *) calloc(2*size_y*local_size_x, sizeof(double));
    uFil_local = (double *) calloc(2*size_y*local_size_x, sizeof(double));
    vFil_local = (double *) calloc(2*size_y*local_size_x, sizeof(double));
    hPhy_local = (double *) calloc(2*size_y*local_size_x, sizeof(double));
    uPhy_local = (double *) calloc(2*size_y*local_size_x, sizeof(double));
    vPhy_local = (double *) calloc(2*size_y*local_size_x, sizeof(double));

}

void dealloc(void) {
  free(hFil);
  free(uFil);
  free(vFil);
  free(hPhy);
  free(uPhy);
  free(vPhy);
}
