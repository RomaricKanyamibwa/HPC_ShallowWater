#include <stdlib.h>
#include <shalw.h>
#include <parse_args.hpp>
#include <memory.h>
#include <init.h>
#include <forward.h>
#include <export.h>
#include <time.h>   /* chronometrage */
#include <sys/time.h>
#include <mpi.h>

double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
int size_x, size_y, nb_steps;
double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
bool file_export;
std::string export_path;
int my_rank,NP/*Nombre de processeur*/;

int main(int argc, char **argv) {

  parse_args(argc, argv);
  printf("Command line options parsed\n");

  alloc();
  printf("Memory allocated\n");

  gauss_init();
  printf("State initialised\n");
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&NP);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  forward();
  printf("State computed\n");

  dealloc();
  printf("Memory freed\n");
  MPI_Finalize();

  return EXIT_SUCCESS;
}
