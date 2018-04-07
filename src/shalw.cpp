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
int NbCol, NbLi, ligne_colonne,NP_temp,local_size_y;//par bloc



double my_gettimeofday(){
  struct timeval tmp_time;
  gettimeofday(&tmp_time, NULL);
  return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}


int main_bloc(int argc, char **argv)  {

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

  /*NP doit etre une puissance de 2*/
  NP_temp = NP;
  local_size_x = size_x;
  local_size_y = size_y;
  NbCol = 0;
  NbLi = 0;
  ligne_colonne = 0; //0 pour ligne 1 pour colonne
  while(NP_temp!=1)
  {
    NP_temp = NP_temp/2;
    if(ligne_colonne==0) //ligne
    {
      local_size_x = local_size_x/2;
      NbLi++;
    }
    else
    {
      local_size_y = local_size_y/2;
      NbCol++;
    }
    ligne_colonne = (ligne_colonne+1)%2;
  }

  if(my_rank>=NbCol)
    local_size_x++;
  if(my_rank<(NbLi-1)*NbCol)
    local_size_x++;
  if(my_rank%NbCol!=0)
    local_size_y++;
  if((my_rank+1)%NbCol!=0)
    local_size_y++;


  gauss_init();
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

int main(int argc, char **argv) {

  /* Variables liees au chronometrage */
  double debut=0, fin=0;
  root = 0;
  parse_args(argc, argv);
  if(my_rank==0)
    printf("Command line options parsed\n");

  alloc();
  if(my_rank==0)
  {
      printf("Memory allocated\n");
      printf("State initialised\n");
  }
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&NP);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  //if(NP>1)
    local_size_x = (my_rank==0||my_rank==(NP-1))?(size_x/NP+1):(size_x/NP+2);
  //else
    //local_size_x=size_x;
  gauss_init();
  //printf("P#%d:local size x:%d , y:%d\n",my_rank,local_size_x,size_y);
  //printf("P#%d:size x:%d , y:%d\n",my_rank,size_x,size_y);

  if(my_rank==0)
  {
    /* debut du chronometrage */
    debut = my_gettimeofday();
    printf("***************NP:%d***************\n",NP);
  }

  if(size_x%NP!=0)
  {
    fprintf(stderr,"ERROR:NP ne divise pas size_y\n");
    return EXIT_FAILURE;
  }

  forward();

  if(my_rank==0)
  {
     /* fin du chronometrage */
    fin = my_gettimeofday();
    printf("State computed\n");
    printf("#%d-Temps total de calcul : %g seconde(s) \n",my_rank,fin - debut);
  }

  dealloc();
  if(my_rank==0)
    printf("Memory freed\n");

  MPI_Finalize();

  return EXIT_SUCCESS;
}
