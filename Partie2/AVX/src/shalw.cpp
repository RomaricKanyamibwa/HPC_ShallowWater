#include <stdlib.h>
#include <shalw.h>
#include <parse_args.hpp>
#include <memory.h>
#include <init.h>
#include <forward.h>
#include <export.h>
#include <time.h>   /* chronometrage */
#include <sys/time.h>
#include <string.h>
double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
int size_x, size_y, nb_steps;
double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
bool file_export;
std::string export_path;

double my_gettimeofday(){
  struct timeval tmp_time;
  gettimeofday(&tmp_time, NULL);
  return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}


int main(int argc, char **argv) {
  parse_args(argc, argv);
  /* Variables liees au chronometrage */
  double debut=0, fin=0;
  printf("Command line options parsed\n");

  alloc();
  printf("Memory allocated\n");

  gauss_init();
  printf("State initialised\n");
  debut = my_gettimeofday();
  forward();
  fin = my_gettimeofday();
  printf("State computed\n");
  printf("Temps total de calcul : %g seconde(s) \n",fin - debut);
  FILE *perf = fopen("perform.txt", "a+");
  char str[512];
  sprintf(str,"*************** AVX NP:1 ***************\n\
size_x:%d , size_y:%d , nbsteps:%d \n\
Temps total de calcul : %g seconde(s)\n\n",
    size_x,size_y,nb_steps,fin-debut);
  fwrite(str,sizeof(char),strlen(str),perf);

  dealloc();
  printf("Memory freed\n");

  return EXIT_SUCCESS;
}
