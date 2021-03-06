#include <stdlib.h>
#include <shalw.h>
#include <malloc.h>

void alloc(void) {
  if ((hFil = (double *) aligned_alloc(32,2*size_x*size_y*sizeof(double))) == NULL){
        fprintf(stderr, "Error while allocating hfill.\n");
    }
  if ((uFil = (double *) aligned_alloc(32,2*size_x*size_y*sizeof(double))) == NULL){
        fprintf(stderr, "Error while allocating ufill.\n");
    }
  if ((vFil = (double *) aligned_alloc(32,2*size_x*size_y*sizeof(double))) == NULL){
        fprintf(stderr, "Error while allocating vfill.\n");
    }
  if ((hPhy = (double *) aligned_alloc(32,2*size_x*size_y*sizeof(double))) == NULL){
        fprintf(stderr, "Error while allocating hphy.\n");
    }
  if ((vPhy = (double *) aligned_alloc(32,2*size_x*size_y*sizeof(double))) == NULL){
        fprintf(stderr, "Error while allocating vphy.\n");
    }
  if ((uPhy = (double *) aligned_alloc(32,2*size_x*size_y*sizeof(double))) == NULL){
        fprintf(stderr, "Error while allocating uphy.\n");
    }
}

void dealloc(void) {
  free(hFil);
  free(uFil);
  free(vFil);
  free(hPhy);
  free(uPhy);
  free(vPhy);
}
