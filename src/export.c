#include <stdio.h>
#include <shalw.h>

FILE *create_file(void) {
  FILE *f;
  char fname[256];

  sprintf(fname, "%s/shalw_%dx%d_T%d_NP%d.sav", export_path.c_str(), size_x, size_y, nb_steps,NP);
  printf("Fname:%s \n",fname);

  f = fopen(fname, "w+b");

  return f;
}
int i=0;
void export_step(FILE *f, int t) {
	//printf("k1,my_rank %d \n",i++);
	//printf("t=%d \n",t );
  	fwrite((void *)&HFIL(t, 0, 0), sizeof(double), size_x * size_y, f);
	//printf("k2\n");
}

void finalize_export(FILE *f) {
	printf("CLOSE\n");
	if(f==NULL)
		printf("tes trop nul\n");
	if(my_rank==0)
        fclose(f);
  	printf("CLOSE1\n");

}

MPI_File* create_file_mpi(void) {
  MPI_File *f;
  int err = 0;
  char fname[256];

  sprintf(fname, "%s/shalw_%dx%d_T%d_NP%d.sav", export_path.c_str(), size_x, size_y, nb_steps,NP);
  printf("Fname:%s \n",fname);

  err=MPI_File_open(MPI_COMM_WORLD,fname,
                  MPI_MODE_RDWR | MPI_MODE_CREATE
                  ,MPI_INFO_NULL,f);
  if (err)
    {
        MPI_Abort(MPI_COMM_WORLD, 911);
    }
  return f;
}

void export_step_mpi(MPI_File *f, int t) {
    MPI_Status status;
    MPI_Offset offset;
    offset = my_rank*size_x/NP *size_y + (t)*size_x*size_y ;
  	//fwrite((void *)&HFIL(t, 0, 0), sizeof(double), size_x * size_y, f);
  	MPI_File_write_at_all(*f, offset, (void *)&HFIL(t, 0, 0),
                  size_x*size_y, MPI_DOUBLE, &status);
}

void finalize_export_mpi(MPI_File *f) {
	printf("CLOSE\n");
    MPI_File_close(f);
}
