#include <stdlib.h>
#include <shalw.h>

void alloc(void) {
  if(my_rank==0)
  {
    if ((hFil = (double *) aligned_alloc(32,2*size_x*size_y*sizeof(double))) == NULL){
          fprintf(stderr, "Error while allocating hfill.\n");
      }
  }

  if ((hFil_local = (double *) aligned_alloc(32,2*size_x*size_y*sizeof(double))) == NULL){
        fprintf(stderr, "Error while allocating hfill.\n");
    }
  if ((uFil_local = (double *) aligned_alloc(32,2*size_x*size_y*sizeof(double))) == NULL){
        fprintf(stderr, "Error while allocating ufill.\n");
    }
  if ((vFil_local = (double *) aligned_alloc(32,2*size_x*size_y*sizeof(double))) == NULL){
        fprintf(stderr, "Error while allocating vfill.\n");
    }
  if ((hPhy_local = (double *) aligned_alloc(32,2*size_x*size_y*sizeof(double))) == NULL){
        fprintf(stderr, "Error while allocating hphy.\n");
    }
  if ((vPhy_local = (double *) aligned_alloc(32,2*size_x*size_y*sizeof(double))) == NULL){
        fprintf(stderr, "Error while allocating vphy.\n");
    }
  if ((uPhy_local = (double *) aligned_alloc(32,2*size_x*size_y*sizeof(double))) == NULL){
        fprintf(stderr, "Error while allocating uphy.\n");
    }
}
void dealloc(void) {
  free(hFil);
//  free(uFil);
//  free(vFil);
//  free(hPhy);
//  free(uPhy);
//  free(vPhy);
  if(my_rank==0)
  {
    //printf("LOCAL free\n");
    free(hFil_local);
    //printf("LOCAL free2\n");
    free(uFil_local);
    free(vFil_local);
    free(hPhy_local);
    free(uPhy_local);
    free(vPhy_local);
  }
}
