#include <stdio.h>
#include <shalw.h>

FILE *create_file(void) {
  FILE *f;
  char fname[256];

  sprintf(fname, "%s/shalw_%dx%d_T%d.sav", export_path.c_str(), size_x, size_y, nb_steps);

  f = fopen(fname, "w+b");

  return f;
}
int i=0;
void export_step(FILE *f, int t) {
	printf("k1,my_rank %d \n",i++);
	printf("t=%d \n",t );
  	fwrite((void *)&HFIL(t, 0, 0), sizeof(double), size_x * size_y, f);
	printf("k2\n");
}

void finalize_export(FILE *f) {
	printf("CLOSE\n");
	if(f==NULL)
		printf("tes trop nul\n");
	if(my_rank==0)
  fclose(f);
  	printf("CLOSE1\n");

}
