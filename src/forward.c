#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>

#define TAG_MESSAGE
#ifdef TAG_MESSAGE
#define TAG_LAST 0
#define TAG_FIRST 1
#endif /* #ifdef TAG_MESSAGE */

double hFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //HPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return HPHY(t, i, j);
  return HPHY(t - 1, i, j) +
    alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

double uFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //UPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return UPHY(t, i, j);
  return UPHY(t - 1, i, j) +
    alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

double vFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //VPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return VPHY(t, i, j);
  return VPHY(t - 1, i, j) +
    alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
}

double hPhy_forward(int t, int i, int j) {
  double c, d;

  c = 0.;
  if (i > 0)
    c = UPHY(t - 1, i - 1, j);

  d = 0.;
  if (j < size_y - 1)
    d = VPHY(t - 1, i, j + 1);

  return HFIL(t - 1, i, j) -
    dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx +
		 (d - VPHY(t - 1, i, j)) / dy);
}

double uPhy_forward(int t, int i, int j) {
  double b, e, f, g;

  if (i == size_x - 1)
    return 0.;

  b = 0.;
  if (i < size_x - 1)
    b = HPHY(t - 1, i + 1, j);

  e = 0.;
  if (j < size_y - 1)
    e = VPHY(t - 1, i, j + 1);

  f = 0.;
  if (i < size_x - 1)
    f = VPHY(t - 1, i + 1, j);

  g = 0.;
  if (i < size_x - 1 && j < size_y - 1)
    g = VPHY(t - 1, i + 1, j + 1);

  return UFIL(t - 1, i, j) +
    dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) +
	  (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) -
	  (dissip * UFIL(t - 1, i, j)));
}

double vPhy_forward(int t, int i, int j) {
  double c, d, e, f;

  if (j == 0)
    return 0.;

  c = 0.;
  if (j > 0)
    c = HPHY(t - 1, i, j - 1);

  d = 0.;
  if (i > 0 && j > 0)
    d = UPHY(t - 1, i -1, j -1);

  e = 0.;
  if (i > 0)
    e = UPHY(t - 1, i - 1, j);

  f = 0.;
  if (j > 0)
    f = UPHY(t - 1, i, j - 1);

  return VFIL(t - 1, i, j) +
    dt * ((-grav / dy) * (HPHY(t - 1, i, j) - c) -
	  (pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) -
	  (dissip * VFIL(t - 1, i, j)));
}

void forward(void) {
  FILE *file = NULL;
  double svdt = 0.;
  int t = 0,k;
  MPI_Status status;
  int mpi_ret_type;

  if (file_export) {
    file = create_file();
    export_step(file, t);
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

    for(k=-1;k<1;k++)
    {
        if(my_rank!=0)
        {
            mpi_ret_type = MPI_Sendrecv(&HPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST
            ,&HPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&UPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST
            ,&UPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&VPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST
            ,&VPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&UFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST
            ,&UFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&VFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST
            ,&VFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST, MPI_COMM_WORLD,&status);
        }
        if(my_rank!=NP-1)
        {
            mpi_ret_type = MPI_Sendrecv(&HPHY_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_LAST
            ,&HPHY_LOCAL(t + k,local_size_x-1, 0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&UPHY_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_LAST
            ,&UPHY_LOCAL(t + k,local_size_x-1, 0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&VPHY_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_LAST
            ,&VPHY_LOCAL(t + k,local_size_x-1, 0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&UFIL_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_LAST
            ,&UFIL_LOCAL(t + k,local_size_x-1, 0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&VFIL_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_LAST
            ,&VFIL_LOCAL(t + k,local_size_x-1, 0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST, MPI_COMM_WORLD,&status);
        }
    }

    for (int j = 0; j < size_y; j++) {
      for (int i = 0; i < size_x; i++) {
        HPHY_LOCAL(t, i, j) = hPhy_forward(t, i, j);
        UPHY_LOCAL(t, i, j) = uPhy_forward(t, i, j);
        VPHY_LOCAL(t, i, j) = vPhy_forward(t, i, j);
        HFIL_LOCAL(t, i, j) = hFil_forward(t, i, j);
        UFIL_LOCAL(t, i, j) = uFil_forward(t, i, j);
        VFIL_LOCAL(t, i, j) = vFil_forward(t, i, j);
      }
    }

    if (file_export) {
      export_step(file, t);
    }

    if (t == 2) {
      dt = svdt;
    }
  }

  if (file_export) {
    finalize_export(file);
  }
}
