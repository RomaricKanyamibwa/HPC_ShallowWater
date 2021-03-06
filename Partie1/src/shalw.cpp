#include <stdlib.h>
#include <shalw.h>
#include <parse_args.hpp>
#include <memory.h>
#include <init.h>
#include <forward.h>
#include <export.h>
#include <time.h>   /* chronometrage */
#include <sys/time.h>

double *hFil;//, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
double *hFil_local, *uFil_local, *vFil_local, *hPhy_local, *uPhy_local, *vPhy_local;
int size_x, size_y, nb_steps,local_size_x;
double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
bool file_export,decomp_bloc,non_block_comm,pararel_IO,non_block_pararel_IO;
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

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&NP);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

  //parse_args(argc, argv);
  if(my_rank==0)
    printf("Command line options parsed\n");

  alloc();
  if(my_rank==0)
  {
    printf("Memory allocated\n");
    printf("State initialised\n");
  }
  /*NP doit etre une puissance de 2*/
  NP_temp = NP;
  local_size_x = size_x;
  local_size_y = size_y;
  NbCol = 1;
  NbLi = 1;
  ligne_colonne = 0; //0 pour ligne 1 pour colonne
  //printf("Decomposition par bloc\n");
  while(NP_temp!=1)
  {
    NP_temp = NP_temp/2;
    if(ligne_colonne==0) //ligne
    {
      local_size_y = local_size_y/2;
      NbCol*=2;
    }
    else
    {
      local_size_x = local_size_x/2;
      NbLi*=2;
    }
    ligne_colonne = (ligne_colonne+1)%2;
  }
  //printf("End Calcul de size_x et size_y locaux\n");
  if(my_rank>=NbCol)//tous les blocs sauf ceux de la 1ere ligne
    local_size_x++;
  if(my_rank<(NbLi-1)*NbCol)//tous les blocs sauf ceux de la derniere ligne
    local_size_x++;
  if(my_rank%NbCol!=0)//tous les blocs sauf ceux de la colonne la plus a gauche
    local_size_y++;
  if((my_rank+1)%NbCol!=0)//tous les blocs sauf ceux de la colonne la plus a droire
    local_size_y++;


  gauss_init_bloc();
  printf("P#%d:local size x:%d , y:%d\n",my_rank,local_size_x,local_size_y);
  if(my_rank==0)
  {
      printf("P#%d:size x:%d , y:%d\n",my_rank,size_x,size_y);
      printf("P#%d:Nbline:%d , Nbcolonne:%d\n",my_rank,NbLi,NbCol);
  }

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
  forward_bloc();
  printf("State computed\n");

  if(my_rank==0)
  {
     /* fin du chronometrage */
    fin = my_gettimeofday();
    printf("#%d-Temps total de calcul : %g seconde(s) \n",my_rank,fin - debut);
    FILE *perf = fopen("perform.txt", "a+");
    char str[512];
    char tmp[128];
    if(non_block_comm)
        sprintf(tmp,"Decomp_Bloc Non-block Mode");
    else
        sprintf(tmp,"Decomp_Bloc Block-Mode");

    if(non_block_pararel_IO)
        sprintf(tmp,"%s Non-Block-MP_IO",tmp);
    else
    {
        if(pararel_IO)
            sprintf(tmp,"%s MP_IO",tmp);
    }
    sprintf(str,"***************NP:%d - %s***************\n\
size_x:%d , size_y:%d , nbsteps:%d \n\
#%d-Temps total de calcul : %g seconde(s)\n\n"
            ,NP,tmp,size_x,size_y,nb_steps,my_rank,fin-debut);
    fwrite(str,sizeof(char),strlen(str),perf);
  }

  dealloc();
  if(my_rank==0)
    printf("Memory freed\n");

  MPI_Finalize();

  return EXIT_SUCCESS;



}

int main(int argc, char **argv) {

  /* Variables liees au chronometrage */
  double debut=0, fin=0;
  root = 0;
  parse_args(argc, argv);
  if(decomp_bloc)
  {
      return main_bloc(argc,argv);
  }
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&NP);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

  //non_block_comm=true;
  if(my_rank==0)
    printf("Command line options parsed\n");

  alloc();
  if(my_rank==0)
  {
      printf("Memory allocated\n");
      printf("State initialised\n");
  }

  if(size_x%NP!=0)
  {
    fprintf(stderr,"ERROR:NP ne divise pas size_y\n");
    return EXIT_FAILURE;
  }

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
    printf("***************NP:%d***************\n",NP);
    debut = my_gettimeofday();
  }

  if(pararel_IO)
    forward_parallel_io();
  else
    forward();

  if(my_rank==0)
  {
     /* fin du chronometrage */
    fin = my_gettimeofday();
    printf("State computed\n");
    printf("#%d-Temps total de calcul : %g seconde(s) \n",my_rank,fin - debut);
    FILE *perf = fopen("perform.txt", "a+");
    char str[512];
    char tmp[64];
    if(non_block_comm)
        sprintf(tmp,"Non-block Mode");
    else
        sprintf(tmp,"Block-Mode");

    if(non_block_pararel_IO)
        sprintf(tmp,"%s Non-Block-MP_IO",tmp);
    else
    {
        if(pararel_IO)
            sprintf(tmp,"%s MP_IO",tmp);
    }
    sprintf(str,"***************NP:%d - %s***************\n\
size_x:%d , size_y:%d , nbsteps:%d \n\
#%d-Temps total de calcul : %g seconde(s)\n\n"
            ,NP,tmp,size_x,size_y,nb_steps,my_rank,fin-debut);
    fwrite(str,sizeof(char),strlen(str),perf);
  }

  dealloc();
  if(my_rank==0)
    printf("Memory freed\n");

  MPI_Finalize();

  return EXIT_SUCCESS;
}
