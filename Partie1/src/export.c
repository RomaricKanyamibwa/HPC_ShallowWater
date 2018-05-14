#include <stdio.h>
#include <shalw.h>
#include <stdlib.h>
#include <errno.h>


FILE *create_file(void) {
  FILE *f;
  char fname[256];
  //int errnum;
  sprintf(fname, "%s/shalw_%dx%d_T%d_NP%d.sav", export_path.c_str(), size_x, size_y, nb_steps,NP);
  //printf("Fname:%s \n",fname);
  f = fopen(fname, "w+b");
  if (f == NULL) {

      //errnum = errno;
      //fprintf(stderr, "Value of errno: %d\n", errno);
      perror("Error opening file");
      //fprintf(stderr, "Error opening file: %s\n", strerror( errnum ));
   }


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
//	if(f==NULL)
//		printf("tes trop nul\n");
	if(my_rank==0)
        fclose(f);
  	//printf("CLOSE1\n");

}

MPI_File* create_file_mpi(MPI_File *f) {
  //MPI_File *f=NULL;
  int err = 0;
  char fname[256];

  sprintf(fname, "%s/shalw_%dx%d_T%d_NP%d.sav", export_path.c_str(), size_x, size_y, nb_steps,NP);
  //printf("Fname:%s \n",fname);
  //printf("Begin open\n");
  err=MPI_File_open(MPI_COMM_WORLD,fname,
                  MPI_MODE_RDWR | MPI_MODE_CREATE
                  ,MPI_INFO_NULL,f);
//  printf("End open\n Error message%d\n",err);
  if (err)
    {
        char* string=(char*)calloc(512,sizeof(char));
        int resultlen;
        MPI_Error_string(err,string,&resultlen);
        printf("Process #%d failed to open file %s\nERROR:%s\n",my_rank,fname,string);
        MPI_Abort(MPI_COMM_WORLD, 911);
    }
  return f;
}

void export_step_mpi(MPI_File *f, int t) {
    MPI_Status status;
    MPI_Offset offset;
    if(decomp_bloc)
    {
        offset = my_rank*size_x/NbLi *size_y/NbCol*sizeof(double) + (t)*size_x*size_y*sizeof(double);
        //un MPI_File_set_view a definir
        MPI_File_write_at_all(*f, offset, (void *)&HFIL_LOCAL(t,(my_rank!=0), 0),
                              size_x/NbLi*size_y/NbCol, MPI_DOUBLE, &status);
        return;
    }
    else
    {
        offset = my_rank*size_x/NP *size_y*sizeof(double) + (t)*size_x*size_y*sizeof(double);
    }
  	//fwrite((void *)&HFIL(t, 0, 0), sizeof(double), size_x * size_y, f);
  	MPI_File_write_at_all(*f, offset, (void *)&HFIL_LOCAL(t,(my_rank!=0), 0),
                  size_x/NP*size_y, MPI_DOUBLE, &status);
}


void export_step_mpi_begin(MPI_File *f, int t) {
    MPI_Offset offset;
    if(decomp_bloc)
    {
        //un MPI_File_set_view a defini
        offset = my_rank*size_x/NbLi *size_y/NbCol*sizeof(double) + (t)*size_x*size_y*sizeof(double);
        MPI_File_write_at_all_begin(*f, offset, (void *)&HFIL_LOCAL(t,(my_rank!=0), 0),
                  size_x/NbLi*size_y/NbCol, MPI_DOUBLE);
        return;
    }
    else
    {
        offset = my_rank*size_x/NP *size_y*sizeof(double) + (t)*size_x*size_y*sizeof(double);
    }
  	//fwrite((void *)&HFIL(t, 0, 0), sizeof(double), size_x * size_y, f);
  	MPI_File_write_at_all_begin(*f, offset, (void *)&HFIL_LOCAL(t,(my_rank!=0), 0),
                  size_x/NP*size_y, MPI_DOUBLE);
}


void export_step_mpi_end(MPI_File *f, int t) {
    MPI_Status status;
  	//fwrite((void *)&HFIL(t, 0, 0), sizeof(double), size_x * size_y, f);
  	MPI_File_write_at_all_end(*f, (void *)&HFIL_LOCAL(t,(my_rank!=0), 0),&status);
}

void finalize_export_mpi(MPI_File *f) {
	printf("CLOSE\n");
    MPI_File_close(f);
}
