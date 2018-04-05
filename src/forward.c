#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>

#define TAG_MESSAGE
#ifdef TAG_MESSAGE
#define TAG_FIRST_U_P 0
#define TAG_LAST_U_P 0
#define TAG_FIRST_V_P 0//101
#define TAG_LAST_V_P 0//100
#define TAG_FIRST_H_P 0//201
#define TAG_LAST_H_P 0//200
#define TAG_FIRST_U_F 0//301
#define TAG_LAST_U_F 0//300
#define TAG_FIRST_V_F 0//401
#define TAG_LAST_V_F 0//400
#define TAG_FIRST_H_F 0//501
#define TAG_LAST_H_F 0//500
#endif /* #ifdef TAG_MESSAGE */

double hFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //HPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return HPHY_LOCAL(t, i, j);
  return HPHY_LOCAL(t - 1, i, j) +
    alpha * (HFIL_LOCAL(t - 1, i, j) - 2 * HPHY_LOCAL(t - 1, i, j) + HPHY_LOCAL(t, i, j));
}

double uFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //UPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return UPHY_LOCAL(t, i, j);
  return UPHY_LOCAL(t - 1, i, j) +
    alpha * (UFIL_LOCAL(t - 1, i, j) - 2 * UPHY_LOCAL(t - 1, i, j) + UPHY_LOCAL(t, i, j));
}

double vFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //VPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return VPHY_LOCAL(t, i, j);
  return VPHY_LOCAL(t - 1, i, j) +
    alpha * (VFIL_LOCAL(t - 1, i, j) - 2 * VPHY_LOCAL(t - 1, i, j) + VPHY_LOCAL(t, i, j));
}

double hPhy_forward(int t, int i, int j) {
  double c, d;

  c = 0.;
  if (i > 0)
    c = UPHY_LOCAL(t - 1, i - 1, j);

  d = 0.;
  if (j < size_y - 1)
    d = VPHY_LOCAL(t - 1, i, j + 1);

  return HFIL_LOCAL(t - 1, i, j) -
    dt * hmoy * ((UPHY_LOCAL(t - 1, i, j) - c) / dx +
		 (d - VPHY_LOCAL(t - 1, i, j)) / dy);
}

double uPhy_forward(int t, int i, int j) {
  double b, e, f, g;

  if (i == size_x/NP - 1)
    return 0.;

  b = 0.;
  if (i < size_x/NP - 1)
    b = HPHY_LOCAL(t - 1, i + 1, j);

  e = 0.;
  if (j < size_y - 1)
    e = VPHY_LOCAL(t - 1, i, j + 1);

  f = 0.;
  if (i < size_x/NP - 1)
    f = VPHY_LOCAL(t - 1, i + 1, j);

  g = 0.;
  if (i < size_x/NP - 1 && j < size_y - 1)
    g = VPHY_LOCAL(t - 1, i + 1, j + 1);

  return UFIL_LOCAL(t - 1, i, j) +
    dt * ((-grav / dx) * (b - HPHY_LOCAL(t - 1, i, j)) +
	  (pcor / 4.) * (VPHY_LOCAL(t - 1, i, j) + e + f + g) -
	  (dissip * UFIL_LOCAL(t - 1, i, j)));
}

double vPhy_forward(int t, int i, int j) {
  double c, d, e, f;

  if (j == 0)
    return 0.;

  c = 0.;
  if (j > 0)
    c = HPHY_LOCAL(t - 1, i, j - 1);

  d = 0.;
  if (i > 0 && j > 0)
    d = UPHY_LOCAL(t - 1, i -1, j -1);

  e = 0.;
  if (i > 0)
    e = UPHY_LOCAL(t - 1, i - 1, j);

  f = 0.;
  if (j > 0)
    f = UPHY_LOCAL(t - 1, i, j - 1);

  return VFIL_LOCAL(t - 1, i, j) +
    dt * ((-grav / dy) * (HPHY_LOCAL(t - 1, i, j) - c) -
	  (pcor / 4.) * (d + e + f + UPHY_LOCAL(t - 1, i, j)) -
	  (dissip * VFIL_LOCAL(t - 1, i, j)));
}

void forward(void) {
  FILE *file = NULL;
  double svdt = 0.;
  int t = 0,k;
  MPI_Status status;
  int mpi_ret_type;

  if(my_rank==0)
  {
	  if (file_export) {
	    printf("1\n");
	    file = create_file();
	    printf("2\n");
	    export_step(file, t);
	    printf("3\n");
	  }
  }

  for (t = 1; t < nb_steps; t++) {
    if (t == 1) {
      svdt = dt;
      dt = 0;
    }
    if (t == 2){
      dt = svdt / 2.;
    }
    //double  *uFil_local, *vFil_local, *hPhy_local, *uPhy_local, *vPhy_local;

    for(k=0;k<2;k++)
    {
        if(my_rank!=0)
        {
            mpi_ret_type = MPI_Sendrecv(&HPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_H_P
            ,&HPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&UPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_U_P
            ,&UPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_U_P, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&VPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_V_P
            ,&VPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_V_P, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&UFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_U_F
            ,&UFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_U_F, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&VFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_V_F
            ,&VFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_V_F, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&HFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_H_F
            ,&HFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_F, MPI_COMM_WORLD,&status);

            //printf("P#%d:mpirettype_1%d\n",my_rank, mpi_ret_type);

        }
        if(my_rank!=NP-1)
        {
            mpi_ret_type = MPI_Sendrecv(&HPHY_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_H_P
            ,&HPHY_LOCAL(t + k,local_size_x-1, 0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_H_P, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&UPHY_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_U_P
            ,&UPHY_LOCAL(t + k,local_size_x-1, 0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_U_P, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&VPHY_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_V_P
            ,&VPHY_LOCAL(t + k,local_size_x-1, 0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_V_P, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&UFIL_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_U_F
            ,&UFIL_LOCAL(t + k,local_size_x-1, 0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_U_F, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&VFIL_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_V_F
            ,&VFIL_LOCAL(t + k,local_size_x-1, 0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_V_F, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&HFIL_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_H_F
            ,&HFIL_LOCAL(t + k,local_size_x-1, 0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_H_F, MPI_COMM_WORLD,&status);

            //printf("P#%d:mpirettype_2%d\n",my_rank, mpi_ret_type);
        }
    }

    mpi_ret_type++;
    //printf("P#%d:line%d\n",my_rank,189);
    for (int j = 0; j < size_y; j++) {
      for (int i = 0; i < size_x/NP; i++) {
          if(my_rank==0)
          {
            HPHY_LOCAL(t, i, j) = hPhy_forward(t, i, j);
            UPHY_LOCAL(t, i, j) = uPhy_forward(t, i, j);
            VPHY_LOCAL(t, i, j) = vPhy_forward(t, i, j);
            HFIL_LOCAL(t, i, j) = hFil_forward(t, i, j);
            UFIL_LOCAL(t, i, j) = uFil_forward(t, i, j);
            VFIL_LOCAL(t, i, j) = vFil_forward(t, i, j);
          }
          else
          {
            HPHY_LOCAL(t, i+1, j) = hPhy_forward(t, i+1, j);
            UPHY_LOCAL(t, i+1, j) = uPhy_forward(t, i+1, j);
            VPHY_LOCAL(t, i+1, j) = vPhy_forward(t, i+1, j);
            HFIL_LOCAL(t, i+1, j) = hFil_forward(t, i+1, j);
            UFIL_LOCAL(t, i+1, j) = uFil_forward(t, i+1, j);
            VFIL_LOCAL(t, i+1, j) = vFil_forward(t, i+1, j);
          }
      }
    }
    printf("---------------------------- Magic The Gathering ----------------------------\n");
    MPI_Gather(&HFIL_LOCAL(t,(my_rank!=0), 0)/*+size_y*(my_rank!=0)*/,size_y*size_x/NP/*(local_size_x-1-1*(my_rank!=0 && my_rank!=NP-1))*/
    ,MPI_DOUBLE,&HFIL(t, 0, 0),size_y*size_x/NP/*(local_size_x-1-1*(my_rank!=0 && my_rank!=NP-1))*/,MPI_DOUBLE,0,MPI_COMM_WORLD);

//    MPI_Gather(&UFIL_LOCAL(t,(my_rank!=0), 0),size_y*size_x/NP
//    ,MPI_DOUBLE,&UFIL(t, 0, 0),size_y*size_x/NP,MPI_DOUBLE,0,MPI_COMM_WORLD);
//
//    MPI_Gather(&VFIL_LOCAL(t,(my_rank!=0), 0),size_y*size_x/NP
//    ,MPI_DOUBLE,&VFIL(t, 0, 0),size_y*size_x/NP,MPI_DOUBLE,0,MPI_COMM_WORLD);
//
//    MPI_Gather(&HPHY_LOCAL(t,(my_rank!=0), 0),size_y*size_x/NP
//    ,MPI_DOUBLE,&HPHY(t, 0, 0),size_y*size_x/NP,MPI_DOUBLE,0,MPI_COMM_WORLD);
//
//    MPI_Gather(&UPHY_LOCAL(t,(my_rank!=0), 0),size_y*size_x/NP
//    ,MPI_DOUBLE,&UPHY(t, 0, 0),size_y*size_x/NP,MPI_DOUBLE,0,MPI_COMM_WORLD);
//
//    MPI_Gather(&VPHY_LOCAL(t,(my_rank!=0), 0),size_y*size_x/NP
//    ,MPI_DOUBLE,&VPHY(t, 0, 0),size_y*size_x/NP,MPI_DOUBLE,0,MPI_COMM_WORLD);
    printf("---------------------------- End of The Gathering ----------------------------\n");
//    printf("P#%d:---memcpy\n",my_rank);
//    memcpy(&HFIL(t,0,0),&HFIL_LOCAL(t,(my_rank!=0),0),size_x*size_y);
//    memcpy(&VFIL(t,0,0),&VFIL_LOCAL(t,(my_rank!=0),0),size_x*size_y);
//    memcpy(&UFIL(t,0,0),&UFIL_LOCAL(t,(my_rank!=0),0),size_x*size_y);
//    memcpy(&HPHY(t,0,0),&HPHY_LOCAL(t,(my_rank!=0),0),size_x*size_y);
//    memcpy(&HFIL(t,0,0),&UPHY_LOCAL(t,(my_rank!=0),0),size_x*size_y);
//    memcpy(&HFIL(t,0,0),&HFIL_LOCAL(t,(my_rank!=0),0),size_x*size_y);
//    hFil[(0)+(0)*size_y +((t)%2) * size_x * size_y ] = hFil_local[(0)+(0)*size_y +((t)%2) *size_x * size_y];
//    uFil[(0)+(0)*size_y +((t)%2) * size_x * size_y ] = uFil_local[(0)+(0)*size_y +((t)%2) *size_x * size_y];
//    vFil[(0)+(0)*size_y +((t)%2) * size_x * size_y ] = vFil_local[(0)+(0)*size_y +((t)%2) *size_x * size_y];
//    uPhy[(0)+(0)*size_y +((t)%2) * size_x * size_y ] = uPhy_local[(0)+(0)*size_y +((t)%2) *size_x * size_y];
//    vPhy[(0)+(0)*size_y +((t)%2) * size_x * size_y ] = vPhy_local[(0)+(0)*size_y +((t)%2) *size_x * size_y];
//    hPhy[(0)+(0)*size_y +((t)%2) * size_x * size_y ] = hPhy_local[(0)+(0)*size_y +((t)%2) *size_x * size_y];
	if(my_rank==0)
	{
	    if (file_export) {
	      export_step(file, t);
	    }
	    printf("export_step\n");
	}
    if (t == 2) {
      dt = svdt;
    }
  }
  if(my_rank==0)
  if (file_export) {
  	printf("finalize_export\n");
    finalize_export(file);
  }
}
