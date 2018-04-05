#include <stdlib.h>
#include <shalw.h>
#include <parse_args.hpp>
#include <memory.h>
#include <init.h>
#include <forward.h>
#include <export.h>
#include <time.h>   /* chronometrage */
#include <sys/time.h>

double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
double *hFil_local, *uFil_local, *vFil_local, *hPhy_local, *uPhy_local, *vPhy_local;
int size_x, size_y, nb_steps,local_size_x;
double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
bool file_export;
std::string export_path;
int my_rank,NP/*Nombre de processeur*/,root;



double my_gettimeofday(){
  struct timeval tmp_time;
  gettimeofday(&tmp_time, NULL);
  return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}


int main(int argc, char **argv) {

  /* Variables liees au chronometrage */
  double debut=0, fin=0;
  root = 0;
  parse_args(argc, argv);
  printf("Command line options parsed\n");

  alloc();
  printf("Memory allocated\n");
  printf("State initialised\n");
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&NP);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  gauss_init();

  //if(NP>1)
    local_size_x = (my_rank==0||my_rank==(NP-1))?(size_x/NP+1):(size_x/NP+2);
  //else
    //local_size_x=size_x;
  printf("P#%d:local size x:%d , y:%d\n",my_rank,local_size_x,size_y);
  printf("P#%d:size x:%d , y:%d\n",my_rank,size_x,size_y);

  if(my_rank==0)
  {
    /* debut du chronometrage */
    debut = my_gettimeofday();
    printf("NP:%d\n",NP);
  }

  if(size_x%NP!=0)
  {
    fprintf(stderr,"ERROR:NP ne divise pas size_y\n");
    return EXIT_FAILURE;
  }

  forward();
  printf("State computed\n");

  if(my_rank==0)
  {
     /* fin du chronometrage */
    fin = my_gettimeofday();
    printf("#%d-Temps total de calcul : %g seconde(s) \n",my_rank,fin - debut);
  }

  dealloc();
  printf("Memory freed\n");

  MPI_Finalize();

  return EXIT_SUCCESS;
}
