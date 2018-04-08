#include <stdio.h>

FILE *create_file(void);
void export_step(FILE*, int);
void finalize_export(FILE*);

MPI_File *create_file_mpi(MPI_File*);
void export_step_mpi(MPI_File*, int);
void finalize_export_mpi(MPI_File*);
void export_step_mpi_begin(MPI_File *f, int t);
void export_step_mpi_end(MPI_File *f, int t);
