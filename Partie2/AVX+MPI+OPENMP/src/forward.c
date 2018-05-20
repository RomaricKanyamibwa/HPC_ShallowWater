#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>
#include <stdlib.h>
#include <string.h>

#define TAG_MESSAGE
#ifdef TAG_MESSAGE
#define TAG_FIRST_U_P 1
#define TAG_LAST_U_P 0
#define TAG_FIRST_V_P 101
#define TAG_LAST_V_P 100
#define TAG_FIRST_H_P 201
#define TAG_LAST_H_P 200
#define TAG_FIRST_U_F 301
#define TAG_LAST_U_F 300
#define TAG_FIRST_V_F 401
#define TAG_LAST_V_F 400
#define TAG_FIRST_H_F 501
#define TAG_LAST_H_F 500
/**bande (horizonale) haut et bas**/
#define TAG_BLOC_HOR_FIRST_U_P 1
#define TAG_BLOC_HOR_LAST_U_P 0
#define TAG_BLOC_HOR_FIRST_V_P 101
#define TAG_BLOC_HOR_LAST_V_P 100
#define TAG_BLOC_HOR_FIRST_H_P 201
#define TAG_BLOC_HOR_LAST_H_P 200
#define TAG_BLOC_HOR_FIRST_U_F 301
#define TAG_BLOC_HOR_LAST_U_F 300
#define TAG_BLOC_HOR_FIRST_V_F 401
#define TAG_BLOC_HOR_LAST_V_F 400
#define TAG_BLOC_HOR_FIRST_H_F 501
#define TAG_BLOC_HOR_LAST_H_F 500

/**bande (verticale) du gauche et droite**/
#define TAG_BLOC_VER_FIRST_U_P 601
#define TAG_BLOC_VER_LAST_U_P 600
#define TAG_BLOC_VER_FIRST_V_P 701
#define TAG_BLOC_VER_LAST_V_P 700
#define TAG_BLOC_VER_FIRST_H_P 801
#define TAG_BLOC_VER_LAST_H_P 800
#define TAG_BLOC_VER_FIRST_U_F 901
#define TAG_BLOC_VER_LAST_U_F 900
#define TAG_BLOC_VER_FIRST_V_F 1001
#define TAG_BLOC_VER_LAST_V_F 1000
#define TAG_BLOC_VER_FIRST_H_F 1101
#define TAG_BLOC_VER_LAST_H_F 1100

#define TAG_GATHER 55
#define TAG_REQ 56

#endif /* #ifdef TAG_MESSAGE */



void hFil_forward_vect(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //HPHY_LOCAL(t - 1, i, j) est encore nul

  __m256d hPhy_vect = _mm256_load_pd(&HPHY_LOCAL(t, i, j*4));

  if (t <= 2)
  {
    _mm256_store_pd(&HFIL_LOCAL(t, i, j*4),hPhy_vect);//return HPHY_LOCAL(t, i, j);
    return;
  }

  const __m256d alpha_vect = _mm256_set1_pd(alpha);
  const __m256d s2 = _mm256_set1_pd(-2);

  __m256d hFil_vect = _mm256_load_pd(&HFIL_LOCAL(t - 1, i, j*4));
  __m256d hPhy_vect2 = _mm256_load_pd(&HPHY_LOCAL(t - 1, i, j*4));
  __m256d res =  _mm256_mul_pd(hPhy_vect2, s2);

  res = _mm256_add_pd(res, hFil_vect);
  res = _mm256_add_pd(res, hPhy_vect);
  res = _mm256_mul_pd(res, alpha_vect);
  res = _mm256_add_pd(res, hPhy_vect2);
  _mm256_store_pd(&HFIL_LOCAL(t, i, j*4),res);
  return;
//  HPHY_LOCAL(t - 1, i, j) +
//    alpha * (HFIL_LOCAL(t - 1, i, j) - 2 * HPHY_LOCAL(t - 1, i, j) + HPHY_LOCAL(t, i, j));
}

void uFil_forward_vect(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //UPHY_LOCAL(t - 1, i, j) est encore nul
  __m256d uPhy_vect = _mm256_load_pd(&UPHY_LOCAL(t, i, j*4));

  if (t <= 2)
  {
    _mm256_store_pd(&UFIL_LOCAL(t, i, j*4),uPhy_vect);
  }
  const __m256d alpha_vect = _mm256_set1_pd(alpha);
  const __m256d s2 = _mm256_set1_pd(-2);

  __m256d uPhy_vect2 = _mm256_load_pd(&UPHY_LOCAL(t-1, i, j*4));
  __m256d uFil_vect  = _mm256_load_pd(&UFIL_LOCAL(t-1, i, j*4));

  __m256d res = _mm256_mul_pd(uPhy_vect2, s2);
  res = _mm256_add_pd(res, uFil_vect);
  res = _mm256_add_pd(res, uPhy_vect);
  res = _mm256_mul_pd(res, alpha_vect);
  res = _mm256_add_pd(res, uPhy_vect2);

  _mm256_store_pd(&UFIL_LOCAL(t, i, j*4),res);
  //UPHY_LOCAL(t - 1, i, j) +alpha * (UFIL_LOCAL(t - 1, i, j) - 2 * UPHY_LOCAL(t - 1, i, j) + UPHY_LOCAL(t, i, j))
  return;
}


void vFil_forward_vect(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //VPHY_LOCAL(t - 1, i, j) est encore nul
  __m256d vPhy_vect = _mm256_load_pd(&VPHY_LOCAL(t, i, j*4));

  if (t <= 2)
  {
    _mm256_store_pd(&VFIL_LOCAL(t, i, j*4),vPhy_vect);
    return;//VPHY_LOCAL(t, i, j);
  }

  const __m256d alpha_vect = _mm256_set1_pd(alpha);
  const __m256d s2 = _mm256_set1_pd(-2);

  __m256d vPhy_vect2 = _mm256_load_pd(&VPHY_LOCAL(t-1, i, j*4));
  __m256d vFil_vect  = _mm256_load_pd(&VFIL_LOCAL(t-1, i, j*4));

  __m256d res = _mm256_mul_pd(vPhy_vect2, s2);
  res = _mm256_add_pd(res, vFil_vect);
  res = _mm256_add_pd(res, vPhy_vect);
  res = _mm256_mul_pd(res, alpha_vect);
  res = _mm256_add_pd(res, vPhy_vect2);

  _mm256_store_pd(&VFIL_LOCAL(t, i, j*4),res);
  //VPHY_LOCAL(t - 1, i, j) + alpha * (VFIL_LOCAL(t - 1, i, j) - 2 * VPHY_LOCAL(t - 1, i, j) + VPHY_LOCAL(t, i, j))
  return;
}


void hPhy_forward_vect(int t, int i, int j) {
  __m256d c, d;

  if (i > 0)
    c = _mm256_load_pd(&UPHY_LOCAL(t - 1, i - 1, j*4));
  else
    c=_mm256_setzero_pd();//_mm256_set1_pd(0.0);

  if (j < size_y/4 - 1)
    d = _mm256_loadu_pd(&VPHY_LOCAL(t - 1, i, j*4 + 1));
  else
    d = _mm256_set_pd(VPHY_LOCAL(t - 1, i, j*4 + 1),VPHY_LOCAL(t - 1, i, j*4 + 2),VPHY_LOCAL(t - 1, i, j*4 + 3),0.0);

  __m256d hFil_vect = _mm256_load_pd(&HFIL_LOCAL(t-1, i, j*4));
  __m256d uPhy_vect = _mm256_load_pd(&UPHY_LOCAL(t-1, i, j*4));
  __m256d vPhy_vect = _mm256_load_pd(&VPHY_LOCAL(t-1, i, j*4));

  __m256d res=_mm256_sub_pd(uPhy_vect,c);//(UPHY_LOCAL(t - 1, i, j) - c
  __m256d res2=_mm256_sub_pd(d,vPhy_vect);//d - VPHY_LOCAL(t - 1, i, j)
  res2 = _mm256_div_pd(res2,_mm256_set1_pd(dy*1.0));//(d - VPHY_LOCAL(t - 1, i, j)) / dy)
  res = _mm256_div_pd(res,_mm256_set1_pd(dx*1.0));//(UPHY_LOCAL(t - 1, i, j) - c) / dx
  res = _mm256_add_pd(res,res2);//(UPHY_LOCAL(t - 1, i, j) - c) / dx +(d - VPHY_LOCAL(t - 1, i, j)) / dy
  res = _mm256_mul_pd(_mm256_set1_pd(-dt*hmoy*1.0),res);//-dt * hmoy * ((UPHY_LOCAL(t - 1, i, j) - c) / dx +(d - VPHY_LOCAL(t - 1, i, j)) / dy)
  res = _mm256_add_pd(hFil_vect,res);//HFIL_LOCAL(t - 1, i, j) -dt * hmoy * ((UPHY_LOCAL(t - 1, i, j) - c) / dx +(d - VPHY_LOCAL(t - 1, i, j)) / dy)
  _mm256_store_pd(&HPHY_LOCAL(t, i, j*4),res);

  return;
}


void uPhy_forward_vect(int t, int i, int j) {
  __m256d b, e, f, g;

  if (i == size_x - 1)
  {
    _mm256_store_pd(&UPHY_LOCAL(t, i, j*4),_mm256_setzero_pd());

    return;
  }

  if (i < size_x - 1)
    b = _mm256_load_pd(&HPHY_LOCAL(t - 1, i + 1, j*4));
  else
    b=_mm256_setzero_pd();

  if (j < size_y/4 - 1)
    e = _mm256_loadu_pd(&VPHY_LOCAL(t - 1, i, j*4 + 1));
  else
    e=_mm256_set_pd(VPHY_LOCAL(t - 1, i, j*4 + 1),VPHY_LOCAL(t - 1, i, j*4 + 2),VPHY_LOCAL(t - 1, i, j*4 + 3),0.0);

  if (i < size_x - 1)
    f = _mm256_load_pd(&VPHY_LOCAL(t - 1, i + 1, j*4));
  else
    f = _mm256_setzero_pd();

  if (i < size_x - 1 && j < size_y/4 - 1)
    g = _mm256_loadu_pd(&VPHY_LOCAL(t - 1, i + 1, j*4 + 1));
  else
  {
    if(i < size_x -1)
        g = _mm256_set_pd(VPHY_LOCAL(t - 1, i + 1, j*4 + 1),
        VPHY_LOCAL(t - 1, i + 1, j*4 + 2),
        VPHY_LOCAL(t - 1, i + 1, j*4 + 3),0.0);
    else
        g = _mm256_setzero_pd();
  }
    //_mm256_set_pd(VPHY_LOCAL(t - 1, i + 1, j*4 + 1),VPHY_LOCAL(t - 1, i + 1, j*4 + 2),VPHY_LOCAL(t - 1, i + 1, j*4 + 3),0.0);

  __m256d uFil_vect=_mm256_load_pd(&UFIL_LOCAL(t - 1, i, j*4));
  __m256d hPhy_vect=_mm256_load_pd(&HPHY_LOCAL(t - 1, i, j*4));
  __m256d vPhy_vect=_mm256_load_pd(&VPHY_LOCAL(t - 1, i, j*4));

  __m256d res=_mm256_sub_pd(b,hPhy_vect);//b - HPHY_LOCAL(t - 1, i, j)
  __m256d res2=_mm256_add_pd(e,_mm256_add_pd(f,g));//e + f + g
  res2=_mm256_mul_pd(_mm256_set1_pd(pcor / 4. ),_mm256_add_pd(vPhy_vect,res2));
  //(pcor / 4.) * (VPHY_LOCAL(t - 1, i, j) + e + f + g)
  res=_mm256_fmadd_pd(_mm256_set1_pd(-grav / dx*1.0),res,res2);
  //(-grav / dx) * (b - HPHY_LOCAL(t - 1, i, j)) + (pcor / 4.) * (VPHY_LOCAL(t - 1, i, j) + e + f + g)
  res=_mm256_sub_pd(res,_mm256_mul_pd(_mm256_set1_pd(dissip),uFil_vect));
  //(-grav / dx) * (b - HPHY_LOCAL(t - 1, i, j)) + (pcor / 4.) * (VPHY_LOCAL(t - 1, i, j) + e + f + g) -
  //(dissip * UFIL_LOCAL(t - 1, i, j))
  res=_mm256_mul_pd(_mm256_set1_pd(dt),res);//dt * ((-grav / dx) * (b - HPHY_LOCAL(t - 1, i, j)) +
  //(pcor / 4.) * (VPHY_LOCAL(t - 1, i, j) + e + f + g) -
  //(dissip * UFIL_LOCAL(t - 1, i, j)))
  res=_mm256_add_pd(res,uFil_vect);
//  UFIL_LOCAL(t - 1, i, j) +
//    dt * ((-grav / dx) * (b - HPHY_LOCAL(t - 1, i, j)) +
//    (pcor / 4.) * (VPHY_LOCAL(t - 1, i, j) + e + f + g) -
//    (dissip * UFIL_LOCAL(t - 1, i, j)))
  _mm256_store_pd(&UPHY_LOCAL(t, i, j*4),res);


  return;
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

void vPhy_forward_vect(int t, int i, int j) {
  __m256d c, d, e, f;

  if (j == 0)
  {

    VPHY_LOCAL(t, i, 0)=0.0;
    VPHY_LOCAL(t, i, 1)=vPhy_forward(t,i,1);
    VPHY_LOCAL(t, i, 2)=vPhy_forward(t,i,2);
    VPHY_LOCAL(t, i, 3)=vPhy_forward(t,i,3);

    return;
  }

  if (j > 0)
    c = _mm256_loadu_pd(&HPHY_LOCAL(t - 1, i, j*4 - 1));
  else
    c = _mm256_set_pd(0.0,HPHY_LOCAL(t - 1, i, 0),HPHY_LOCAL(t - 1, i , 1),HPHY_LOCAL(t - 1, i, 2));

  if (i > 0 && j > 0)
    d = _mm256_loadu_pd(&UPHY_LOCAL(t - 1, i - 1, j*4 -1));
  else
  {
    if(i>0)
        d = _mm256_set_pd(0.0,UPHY_LOCAL(t - 1, i - 1, 0),UPHY_LOCAL(t - 1, i - 1 , 1),UPHY_LOCAL(t - 1, i - 1, 2));
    else
        d = _mm256_setzero_pd();
  }

  if (i > 0)
    e = _mm256_load_pd(&UPHY_LOCAL(t - 1, i - 1, j*4));
  else
    e = _mm256_setzero_pd();

  if (j > 0)
    f = _mm256_loadu_pd(&UPHY_LOCAL(t - 1, i, j*4 - 1));
  else
    f = _mm256_set_pd(0.0,UPHY_LOCAL(t - 1, i, 0),UPHY_LOCAL(t - 1, i , 1),UPHY_LOCAL(t - 1, i, 2));

  __m256d vFil_vect=_mm256_load_pd(&VFIL_LOCAL(t - 1, i, j*4));
  __m256d hPhy_vect=_mm256_load_pd(&HPHY_LOCAL(t - 1, i, j*4));
  __m256d uPhy_vect=_mm256_load_pd(&UPHY_LOCAL(t - 1, i, j*4));

  __m256d res=_mm256_sub_pd(hPhy_vect,c);//HPHY_LOCAL(t - 1, i, j) - c
  __m256d res2=_mm256_add_pd(d,_mm256_add_pd(e,f));//d + e + f
  res2=_mm256_mul_pd(_mm256_set1_pd(pcor / 4. ),_mm256_add_pd(uPhy_vect,res2));
  //(pcor / 4.) * (d + e + f + UPHY_LOCAL(t - 1, i, j))
  res=_mm256_fmsub_pd(_mm256_set1_pd(-grav / dy*1.0),res,res2);
  //(-grav / dy) * (HPHY_LOCAL(t - 1, i, j) - c) - (pcor / 4.) * (d + e + f + UPHY_LOCAL(t - 1, i, j))
  res=_mm256_sub_pd(res,_mm256_mul_pd(_mm256_set1_pd(dissip),vFil_vect));
  //(-grav / dy) * (HPHY_LOCAL(t - 1, i, j) - c) - (pcor / 4.) * (d + e + f + UPHY_LOCAL(t - 1, i, j)) -
  //(dissip * UFIL_LOCAL(t - 1, i, j))
  res=_mm256_mul_pd(_mm256_set1_pd(dt),res);
  // dt * ((-grav / dy) * (HPHY_LOCAL(t - 1, i, j) - c) -
  // (pcor / 4.) * (d + e + f + UPHY_LOCAL(t - 1, i, j)) -
  // (dissip * VFIL_LOCAL(t - 1, i, j)))
  res=_mm256_add_pd(res,vFil_vect);
  // VFIL_LOCAL(t - 1, i, j) +
  // dt * ((-grav / dy) * (HPHY_LOCAL(t - 1, i, j) - c) -
  // (pcor / 4.) * (d + e + f + UPHY_LOCAL(t - 1, i, j)) -
  // (dissip * VFIL_LOCAL(t - 1, i, j)));
  _mm256_store_pd(&VPHY_LOCAL(t, i, j*4),res);

  return;
}

void forward_bloc(void) {
  FILE *file = NULL;
  double svdt = 0.;
  int t = 0,k = 0;
  MPI_Status status;
  MPI_Status stats[48];
  MPI_Request reqs[48];

  /**buffer for send and recv**/
  double* hphy_send=(double *) calloc(size_x/NbLi,sizeof(double));
  double* hphy_recv=(double *) calloc(size_x/NbLi,sizeof(double));

  double* uphy_send=(double *) calloc(size_x/NbLi,sizeof(double));
  double* uphy_recv=(double *) calloc(size_x/NbLi,sizeof(double));

  double* vphy_send=(double *) calloc(size_x/NbLi,sizeof(double));
  double* vphy_recv=(double *) calloc(size_x/NbLi,sizeof(double));

  double* hfil_send=(double *) calloc(size_x/NbLi,sizeof(double));
  double* hfil_recv=(double *) calloc(size_x/NbLi,sizeof(double));

  double* ufil_send=(double *) calloc(size_x/NbLi,sizeof(double));
  double* ufil_recv=(double *) calloc(size_x/NbLi,sizeof(double));

  double* vfil_send=(double *) calloc(size_x/NbLi,sizeof(double));
  double* vfil_recv=(double *) calloc(size_x/NbLi,sizeof(double));

  //unsigned char connect_msg=0;
  double* hfil_buff_send=(double *) calloc(size_x/NbLi*size_y/NbCol,sizeof(double));
  double* hfil_buff_recv=(double *) calloc(size_x*size_y,sizeof(double));

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
  if(my_rank==0)
    printf("P#%d start de decomp par bloc\n",my_rank);
  for (t = 1; t < nb_steps; t++) {

    if (t == 1) {
      svdt = dt;
      dt = 0;
    }
    if (t == 2){
      dt = svdt / 2.;
    }
    #pragma omp parallel
    for (int i = (my_rank>=NbCol); i < local_size_x; i++) {
        for (int j = (my_rank%NbCol!=0); j < local_size_y; j++) {
            hPhy_forward_vect(t,i,j);  
    uPhy_forward_vect(t,i,j);
    vPhy_forward_vect(t,i,j);
    hFil_forward_vect(t,i,j);
    uFil_forward_vect(t,i,j);
    vFil_forward_vect(t,i,j);
      }
    }

    {
        if(non_block_comm)
        {
            /* au dessus en dessous */
            if(my_rank>=NbCol) //on envoie celui du haut sauf ceux sur la premiere ligne
            {
    //            printf("P#%d:Send bande haut\n",my_rank);

                MPI_Isend(&HPHY_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_BLOC_HOR_LAST_H_P,MPI_COMM_WORLD,&reqs[0]);
                MPI_Irecv(&HPHY_LOCAL(t + k,0, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_BLOC_HOR_FIRST_H_P,MPI_COMM_WORLD,&reqs[1]);

                MPI_Isend(&UPHY_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_BLOC_HOR_LAST_U_P,MPI_COMM_WORLD,&reqs[2]);
                MPI_Irecv(&UPHY_LOCAL(t + k,0, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_BLOC_HOR_FIRST_U_P,MPI_COMM_WORLD,&reqs[3]);

                MPI_Isend(&VPHY_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_BLOC_HOR_LAST_V_P,MPI_COMM_WORLD,&reqs[4]);
                MPI_Irecv(&VPHY_LOCAL(t + k,0, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_BLOC_HOR_FIRST_V_P,MPI_COMM_WORLD,&reqs[5]);

                MPI_Isend(&UFIL_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_BLOC_HOR_LAST_U_F,MPI_COMM_WORLD,&reqs[6]);
                MPI_Irecv(&UFIL_LOCAL(t + k,0, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_BLOC_HOR_FIRST_U_F,MPI_COMM_WORLD,&reqs[7]);

                MPI_Isend(&VFIL_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_BLOC_HOR_LAST_V_F,MPI_COMM_WORLD,&reqs[8]);
                MPI_Irecv(&VFIL_LOCAL(t + k,0, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_BLOC_HOR_FIRST_V_F,MPI_COMM_WORLD,&reqs[9]);

                MPI_Isend(&HFIL_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_BLOC_HOR_LAST_H_F,MPI_COMM_WORLD,&reqs[10]);
                MPI_Irecv(&HFIL_LOCAL(t + k,0, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_BLOC_HOR_FIRST_H_F,MPI_COMM_WORLD,&reqs[11]);

                //printf("P#%d:mpirettype_1%d\n",my_rank, mpi_ret_type);
                //printf("P#%d:Final Send bande haut\n",my_rank);

            }
            if(my_rank<NbCol*(NbLi-1)) //envoie celui du bas
            {
                //printf("P#%d:Send bande du bas\n",my_rank);
                MPI_Isend(&HPHY_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE, my_rank+NbCol,TAG_BLOC_HOR_FIRST_H_P,MPI_COMM_WORLD
                             ,&reqs[12-12*(my_rank<NbCol)]);
                MPI_Irecv(&HPHY_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE, my_rank+NbCol,TAG_BLOC_HOR_LAST_H_P,MPI_COMM_WORLD
                             ,&reqs[13-12*(my_rank<NbCol)]);

                MPI_Isend(&UPHY_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE, my_rank+NbCol,TAG_BLOC_HOR_FIRST_U_P,MPI_COMM_WORLD
                             ,&reqs[14-12*(my_rank<NbCol)]);
                MPI_Irecv(&UPHY_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE,my_rank+NbCol,TAG_BLOC_HOR_LAST_U_P,MPI_COMM_WORLD
                             ,&reqs[15-12*(my_rank<NbCol)]);

                MPI_Isend(&VPHY_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE, my_rank+NbCol,TAG_BLOC_HOR_FIRST_V_P,MPI_COMM_WORLD
                             ,&reqs[16-12*(my_rank<NbCol)]);
                MPI_Irecv(&VPHY_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE,my_rank+NbCol,TAG_BLOC_HOR_LAST_V_P,MPI_COMM_WORLD
                             ,&reqs[17-12*(my_rank<NbCol)]);

                MPI_Isend(&UFIL_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE, my_rank+NbCol,TAG_BLOC_HOR_FIRST_U_F,MPI_COMM_WORLD
                             ,&reqs[18-12*(my_rank<NbCol)]);
                MPI_Irecv(&UFIL_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE,my_rank+NbCol,TAG_BLOC_HOR_LAST_U_F,MPI_COMM_WORLD
                             ,&reqs[19-12*(my_rank<NbCol)]);

                MPI_Isend(&VFIL_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE, my_rank+NbCol,TAG_BLOC_HOR_FIRST_V_F,MPI_COMM_WORLD
                             ,&reqs[20-12*(my_rank<NbCol)]);
                MPI_Irecv(&VFIL_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE,my_rank+NbCol,TAG_BLOC_HOR_LAST_V_F,MPI_COMM_WORLD
                             ,&reqs[21-12*(my_rank<NbCol)]);

                MPI_Isend(&HFIL_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE, my_rank+NbCol,TAG_BLOC_HOR_FIRST_H_F,MPI_COMM_WORLD
                             ,&reqs[22-12*(my_rank<NbCol)]);
                MPI_Irecv(&HFIL_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE,my_rank+NbCol,TAG_BLOC_HOR_LAST_H_F,MPI_COMM_WORLD
                             ,&reqs[23-12*(my_rank<NbCol)]);
    //            printf("P#%d:Final Send bande du bas\n",my_rank);
                //printf("P#%d:mpirettype_2%d\n",my_rank, mpi_ret_type);
            }

            /*a droite a gauche */
            if(my_rank%NbCol!=0) //tout les processus sauf ceux qui sont sur la colonne de gauche, on envoie la colonne tout a gauche
            {
    //            printf("P#%d:Send bande gauche\n",my_rank);
                #pragma omp parallel for schedule(static)
                for(int i=0;i<size_x/NbLi;i++){
                    hphy_send[i]=HPHY_LOCAL(t + k,i+(my_rank>=NbCol), 1);
                    uphy_send[i]=UPHY_LOCAL(t + k,i+(my_rank>=NbCol), 1);
                    vphy_send[i]=VPHY_LOCAL(t + k,i+(my_rank>=NbCol), 1);
                    hfil_send[i]=HFIL_LOCAL(t + k,i+(my_rank>=NbCol), 1);
                    ufil_send[i]=UFIL_LOCAL(t + k,i+(my_rank>=NbCol), 1);
                    vfil_send[i]=VFIL_LOCAL(t + k,i+(my_rank>=NbCol), 1);
                }

                MPI_Isend(hphy_send,size_x/NbLi, MPI_DOUBLE, my_rank-1,TAG_BLOC_VER_LAST_H_P,MPI_COMM_WORLD
                ,&reqs[24-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1)))]);
                MPI_Irecv(hphy_recv,size_x/NbLi, MPI_DOUBLE,my_rank-1,TAG_BLOC_VER_FIRST_H_P,MPI_COMM_WORLD
                ,&reqs[25-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1)))]);

                MPI_Isend(uphy_send,size_x/NbLi, MPI_DOUBLE, my_rank-1,TAG_BLOC_VER_LAST_U_P,MPI_COMM_WORLD
                ,&reqs[26-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1)))]);
                MPI_Irecv(uphy_recv,size_x/NbLi, MPI_DOUBLE,my_rank-1,TAG_BLOC_VER_FIRST_U_P,MPI_COMM_WORLD
                ,&reqs[27-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1)))]);


                MPI_Isend(vphy_send,size_x/NbLi, MPI_DOUBLE, my_rank-1,TAG_BLOC_VER_LAST_V_P,MPI_COMM_WORLD
                ,&reqs[28-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1)))]);
                MPI_Irecv(vphy_recv,size_x/NbLi, MPI_DOUBLE,my_rank-1,TAG_BLOC_VER_FIRST_V_P,MPI_COMM_WORLD
                ,&reqs[29-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1)))]);


                MPI_Isend(hfil_send,size_x/NbLi, MPI_DOUBLE, my_rank-1,TAG_BLOC_VER_LAST_H_F,MPI_COMM_WORLD
                ,&reqs[30-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1)))]);
                MPI_Irecv(hfil_recv,size_x/NbLi, MPI_DOUBLE,my_rank-1,TAG_BLOC_VER_FIRST_H_F,MPI_COMM_WORLD
                ,&reqs[31-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1)))]);


                MPI_Isend(ufil_send,size_x/NbLi, MPI_DOUBLE, my_rank-1,TAG_BLOC_VER_LAST_U_F,MPI_COMM_WORLD
                ,&reqs[32-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1)))]);
                MPI_Irecv(ufil_recv,size_x/NbLi, MPI_DOUBLE,my_rank-1,TAG_BLOC_VER_FIRST_U_F,MPI_COMM_WORLD
                ,&reqs[33-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1)))]);


                MPI_Isend(vfil_send,size_x/NbLi, MPI_DOUBLE, my_rank-1,TAG_BLOC_VER_LAST_V_F,MPI_COMM_WORLD
                ,&reqs[34-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1)))]);
                MPI_Irecv(vfil_recv,size_x/NbLi, MPI_DOUBLE,my_rank-1,TAG_BLOC_VER_FIRST_V_F,MPI_COMM_WORLD
                ,&reqs[35-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1)))]);

                //&VFIL_LOCAL(t + k,i+(my_rank>=NbCol), 0)
                #pragma omp parallel for schedule(static)
                for(int i=0;i<size_x/NbLi;i++){
                    HPHY_LOCAL(t + k,i+(my_rank>=NbCol), 0)=hphy_recv[i];
                    UPHY_LOCAL(t + k,i+(my_rank>=NbCol), 0)=uphy_recv[i];
                    VPHY_LOCAL(t + k,i+(my_rank>=NbCol), 0)=vphy_recv[i];
                    HFIL_LOCAL(t + k,i+(my_rank>=NbCol), 0)=hfil_recv[i];
                    UFIL_LOCAL(t + k,i+(my_rank>=NbCol), 0)=ufil_recv[i];
                    VFIL_LOCAL(t + k,i+(my_rank>=NbCol), 0)=vfil_recv[i];
                }

            }
            if((my_rank+1)%NbCol!=0) //tout les processus sauf ceux sur la colonne de droite, on envoiel aligne de droite
            {
    //            printf("P#%d:Send bande droite\n",my_rank);
                #pragma omp parallel for schedule(static)
                for(int i=0;i<size_x/NbLi;i++){
                    hphy_send[i]=HPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2);
                    uphy_send[i]=UPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2);
                    vphy_send[i]=VPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2);
                    hfil_send[i]=HFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2);
                    ufil_send[i]=UFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2);
                    vfil_send[i]=VFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2);
                }

                MPI_Isend(hphy_send,size_x/NbLi, MPI_DOUBLE, my_rank+1,TAG_BLOC_VER_FIRST_H_P,MPI_COMM_WORLD
                ,&reqs[36-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1))+(my_rank%NbCol==0))]);
                MPI_Irecv(hphy_recv,size_x/NbLi, MPI_DOUBLE,my_rank+1,TAG_BLOC_VER_LAST_H_P,MPI_COMM_WORLD
                ,&reqs[37-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1))+(my_rank%NbCol==0))]);

                MPI_Isend(uphy_send,size_x/NbLi, MPI_DOUBLE, my_rank+1,TAG_BLOC_VER_FIRST_U_P,MPI_COMM_WORLD
                ,&reqs[38-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1))+(my_rank%NbCol==0))]);
                MPI_Irecv(uphy_recv,size_x/NbLi, MPI_DOUBLE,my_rank+1,TAG_BLOC_VER_LAST_U_P,MPI_COMM_WORLD
                ,&reqs[39-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1))+(my_rank%NbCol==0))]);

                MPI_Isend(vphy_send,size_x/NbLi, MPI_DOUBLE, my_rank+1,TAG_BLOC_VER_FIRST_V_P,MPI_COMM_WORLD
                ,&reqs[40-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1))+(my_rank%NbCol==0))]);
                MPI_Irecv(vphy_recv,size_x/NbLi, MPI_DOUBLE,my_rank+1,TAG_BLOC_VER_LAST_V_P,MPI_COMM_WORLD
                ,&reqs[41-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1))+(my_rank%NbCol==0))]);

                MPI_Isend(hfil_send,size_x/NbLi, MPI_DOUBLE, my_rank+1,TAG_BLOC_VER_FIRST_H_F,MPI_COMM_WORLD
                ,&reqs[42-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1))+(my_rank%NbCol==0))]);
                MPI_Irecv(hfil_recv,size_x/NbLi, MPI_DOUBLE,my_rank+1,TAG_BLOC_VER_LAST_H_F,MPI_COMM_WORLD
                ,&reqs[43-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1))+(my_rank%NbCol==0))]);

                MPI_Isend(ufil_send,size_x/NbLi, MPI_DOUBLE, my_rank+1,TAG_BLOC_VER_FIRST_U_F,MPI_COMM_WORLD
                ,&reqs[44-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1))+(my_rank%NbCol==0))]);
                MPI_Irecv(ufil_recv,size_x/NbLi, MPI_DOUBLE,my_rank+1,TAG_BLOC_VER_LAST_U_F,MPI_COMM_WORLD
                ,&reqs[45-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1))+(my_rank%NbCol==0))]);

                MPI_Isend(vfil_send,size_x/NbLi, MPI_DOUBLE, my_rank+1,TAG_BLOC_VER_FIRST_V_F,MPI_COMM_WORLD
                ,&reqs[46-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1))+(my_rank%NbCol==0))]);
                MPI_Irecv(vfil_recv,size_x/NbLi, MPI_DOUBLE,my_rank+1,TAG_BLOC_VER_LAST_V_F,MPI_COMM_WORLD
                ,&reqs[47-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1))+(my_rank%NbCol==0))]);

                #pragma omp parallel for schedule(static)
                for(int i=0;i<size_x/NbLi;i++){
                    HPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1)=hphy_recv[i];
                    UPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1)=uphy_recv[i];
                    VPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1)=vphy_recv[i];
                    HFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1)=hfil_recv[i];
                    UFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1)=ufil_recv[i];
                    VFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1)=vfil_recv[i];
                }
            }
        }
        else
        {
            /* au dessus en dessous */
            if(my_rank>=NbCol) //on envoie celui du haut sauf ceux sur la premiere ligne
            {
    //            printf("P#%d:Send bande haut\n",my_rank);

                MPI_Sendrecv(&HPHY_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_BLOC_HOR_LAST_H_P
                ,&HPHY_LOCAL(t + k,0, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_BLOC_HOR_FIRST_H_P, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(&UPHY_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_BLOC_HOR_LAST_U_P
                ,&UPHY_LOCAL(t + k,0, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_BLOC_HOR_FIRST_U_P, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(&VPHY_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_BLOC_HOR_LAST_V_P
                ,&VPHY_LOCAL(t + k,0, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_BLOC_HOR_FIRST_V_P, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(&UFIL_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_BLOC_HOR_LAST_U_F
                ,&UFIL_LOCAL(t + k,0, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_BLOC_HOR_FIRST_U_F, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(&VFIL_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_BLOC_HOR_LAST_V_F
                ,&VFIL_LOCAL(t + k,0, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_BLOC_HOR_FIRST_V_F, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(&HFIL_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_BLOC_HOR_LAST_H_F
                ,&HFIL_LOCAL(t + k,0, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_BLOC_HOR_FIRST_H_F, MPI_COMM_WORLD,&status);

                //printf("P#%d:mpirettype_1%d\n",my_rank, mpi_ret_type);
                //printf("P#%d:Final Send bande haut\n",my_rank);

            }
            if(my_rank<NbCol*(NbLi-1)) //envoie celui du bas
            {
                //printf("P#%d:Send bande du bas\n",my_rank);
                MPI_Sendrecv(&HPHY_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE, my_rank+NbCol,TAG_BLOC_HOR_FIRST_H_P
                             ,&HPHY_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE,my_rank+NbCol,TAG_BLOC_HOR_LAST_H_P, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(&UPHY_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE, my_rank+NbCol,TAG_BLOC_HOR_FIRST_U_P
                             ,&UPHY_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE,my_rank+NbCol,TAG_BLOC_HOR_LAST_U_P, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(&VPHY_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE, my_rank+NbCol,TAG_BLOC_HOR_FIRST_V_P
                             ,&VPHY_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE,my_rank+NbCol,TAG_BLOC_HOR_LAST_V_P, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(&UFIL_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE, my_rank+NbCol,TAG_BLOC_HOR_FIRST_U_F
                             ,&UFIL_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE,my_rank+NbCol,TAG_BLOC_HOR_LAST_U_F, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(&VFIL_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE, my_rank+NbCol,TAG_BLOC_HOR_FIRST_V_F
                             ,&VFIL_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE,my_rank+NbCol,TAG_BLOC_HOR_LAST_V_F, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(&HFIL_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE, my_rank+NbCol,TAG_BLOC_HOR_FIRST_H_F
                             ,&HFIL_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol
                             , MPI_DOUBLE,my_rank+NbCol,TAG_BLOC_HOR_LAST_H_F, MPI_COMM_WORLD,&status);
    //            printf("P#%d:Final Send bande du bas\n",my_rank);
                //printf("P#%d:mpirettype_2%d\n",my_rank, mpi_ret_type);
            }

            /*a droite a gauche */
            if(my_rank%NbCol!=0) //tout les processus sauf ceux qui sont sur la colonne de gauche, on envoie la colonne tout a gauche
            {
    //            printf("P#%d:Send bande gauche\n",my_rank);
                #pragma omp parallel for schedule(static)
                for(int i=0;i<size_x/NbLi;i++){
                    hphy_send[i]=HPHY_LOCAL(t + k,i+(my_rank>=NbCol), 1);
                    uphy_send[i]=UPHY_LOCAL(t + k,i+(my_rank>=NbCol), 1);
                    vphy_send[i]=VPHY_LOCAL(t + k,i+(my_rank>=NbCol), 1);
                    hfil_send[i]=HFIL_LOCAL(t + k,i+(my_rank>=NbCol), 1);
                    ufil_send[i]=UFIL_LOCAL(t + k,i+(my_rank>=NbCol), 1);
                    vfil_send[i]=VFIL_LOCAL(t + k,i+(my_rank>=NbCol), 1);
                }
                MPI_Sendrecv(hphy_send,size_x/NbLi, MPI_DOUBLE, my_rank-1,TAG_BLOC_VER_LAST_H_P
                ,hphy_recv,size_x/NbLi, MPI_DOUBLE,my_rank-1,TAG_BLOC_VER_FIRST_H_P, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(uphy_send,size_x/NbLi, MPI_DOUBLE, my_rank-1,TAG_BLOC_VER_LAST_U_P
                ,uphy_recv,size_x/NbLi, MPI_DOUBLE,my_rank-1,TAG_BLOC_VER_FIRST_U_P, MPI_COMM_WORLD,&status);


                MPI_Sendrecv(vphy_send,size_x/NbLi, MPI_DOUBLE, my_rank-1,TAG_BLOC_VER_LAST_V_P
                ,vphy_recv,size_x/NbLi, MPI_DOUBLE,my_rank-1,TAG_BLOC_VER_FIRST_V_P, MPI_COMM_WORLD,&status);


                MPI_Sendrecv(hfil_send,size_x/NbLi, MPI_DOUBLE, my_rank-1,TAG_BLOC_VER_LAST_H_F
                ,hfil_recv,size_x/NbLi, MPI_DOUBLE,my_rank-1,TAG_BLOC_VER_FIRST_H_F, MPI_COMM_WORLD,&status);


                MPI_Sendrecv(ufil_send,size_x/NbLi, MPI_DOUBLE, my_rank-1,TAG_BLOC_VER_LAST_U_F
                ,ufil_recv,size_x/NbLi, MPI_DOUBLE,my_rank-1,TAG_BLOC_VER_FIRST_U_F, MPI_COMM_WORLD,&status);


                MPI_Sendrecv(vfil_send,size_x/NbLi, MPI_DOUBLE, my_rank-1,TAG_BLOC_VER_LAST_V_F
                ,vfil_recv,size_x/NbLi, MPI_DOUBLE,my_rank-1,TAG_BLOC_VER_FIRST_V_F, MPI_COMM_WORLD,&status);

                //&VFIL_LOCAL(t + k,i+(my_rank>=NbCol), 0)
                #pragma omp parallel for schedule(static)
                for(int i=0;i<size_x/NbLi;i++){
                    HPHY_LOCAL(t + k,i+(my_rank>=NbCol), 0)=hphy_recv[i];
                    UPHY_LOCAL(t + k,i+(my_rank>=NbCol), 0)=uphy_recv[i];
                    VPHY_LOCAL(t + k,i+(my_rank>=NbCol), 0)=vphy_recv[i];
                    HFIL_LOCAL(t + k,i+(my_rank>=NbCol), 0)=hfil_recv[i];
                    UFIL_LOCAL(t + k,i+(my_rank>=NbCol), 0)=ufil_recv[i];
                    VFIL_LOCAL(t + k,i+(my_rank>=NbCol), 0)=vfil_recv[i];
                }

            }
            if((my_rank+1)%NbCol!=0) //tout les processus sauf ceux sur la colonne de droite, on envoiel aligne de droite
            {
    //            printf("P#%d:Send bande droite\n",my_rank);
                #pragma omp parallel for schedule(static)
                for(int i=0;i<size_x/NbLi;i++){
                    hphy_send[i]=HPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2);
                    uphy_send[i]=UPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2);
                    vphy_send[i]=VPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2);
                    hfil_send[i]=HFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2);
                    ufil_send[i]=UFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2);
                    vfil_send[i]=VFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2);
                }
                MPI_Sendrecv(hphy_send,size_x/NbLi, MPI_DOUBLE, my_rank+1,TAG_BLOC_VER_FIRST_H_P
                ,hphy_recv,size_x/NbLi, MPI_DOUBLE,my_rank+1,TAG_BLOC_VER_LAST_H_P, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(uphy_send,size_x/NbLi, MPI_DOUBLE, my_rank+1,TAG_BLOC_VER_FIRST_U_P
                ,uphy_recv,size_x/NbLi, MPI_DOUBLE,my_rank+1,TAG_BLOC_VER_LAST_U_P, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(vphy_send,size_x/NbLi, MPI_DOUBLE, my_rank+1,TAG_BLOC_VER_FIRST_V_P
                ,vphy_recv,size_x/NbLi, MPI_DOUBLE,my_rank+1,TAG_BLOC_VER_LAST_V_P, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(hfil_send,size_x/NbLi, MPI_DOUBLE, my_rank+1,TAG_BLOC_VER_FIRST_H_F
                ,hfil_recv,size_x/NbLi, MPI_DOUBLE,my_rank+1,TAG_BLOC_VER_LAST_H_F, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(ufil_send,size_x/NbLi, MPI_DOUBLE, my_rank+1,TAG_BLOC_VER_FIRST_U_F
                ,ufil_recv,size_x/NbLi, MPI_DOUBLE,my_rank+1,TAG_BLOC_VER_LAST_U_F, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(vfil_send,size_x/NbLi, MPI_DOUBLE, my_rank+1,TAG_BLOC_VER_FIRST_V_F
                ,vfil_recv,size_x/NbLi, MPI_DOUBLE,my_rank+1,TAG_BLOC_VER_LAST_V_F, MPI_COMM_WORLD,&status);

                #pragma omp parallel for schedule(static)
                for(int i=0;i<size_x/NbLi;i++){
                    HPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1)=hphy_recv[i];
                    UPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1)=uphy_recv[i];
                    VPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1)=vphy_recv[i];
                    HFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1)=hfil_recv[i];
                    UFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1)=ufil_recv[i];
                    VFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1)=vfil_recv[i];
                }
            }

        }
    }

    //printf("P#%d:line%d\n",my_rank,189);

    //printf("P#%d:---------------------------- Magic The Gathering ----------------------------\n",my_rank);
    #pragma omp parallel for schedule(static)
    for(int i=0;i<size_x/NbLi;i++)//construction de buffer ligne par ligne
    {
        memcpy(hfil_buff_send+i*size_y/NbCol,&HFIL_LOCAL(t,i+(my_rank>=NbCol), (my_rank%NbCol!=0)),size_y/NbCol*sizeof(double));
    }
    //printf("P#%d-------End for-------\n",my_rank);
    MPI_Gather(hfil_buff_send,size_y/NbCol*size_x/NbLi,MPI_DOUBLE
               ,hfil_buff_recv//&HFIL(t, (my_rank%NbCol)*size_y/NbCol,(my_rank%NbLi)*size_x/NbLi)
               ,size_y/NbCol*size_x/NbLi,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //printf("P#%d-------End gather-------\n",my_rank);
    if(my_rank==0)
    {
        //printf("P#%d:---------------------------- Magic The Gathering ----------------------------\n",my_rank);
        #pragma omp parallel for schedule(static)
        for(int i=0;i<size_x/NbLi;i++)
        {
            for(int j=0;j<NP;j++)
            {
                memcpy(&HFIL(t, i+(j/NbCol)*size_x/NbLi,(j%NbCol)*size_y/NbCol)
                       ,hfil_buff_recv+j*size_y/NbCol*size_x/NbLi+i*size_y/NbCol,size_y/NbCol*sizeof(double));
            }
        }
        //printf("P#%d:---------------------------- End of The Gathering ----------------------------\n",my_rank);
    }
    //printf("P#%d-------End if-------\n",my_rank);
    //printf("P#%d:---------------------------- End of The Gathering ----------------------------\n",my_rank);
   if(my_rank==0)
	{
	    if (file_export) {
	      export_step(file, t);
	    }
	    //printf("export_step\n");
	}
    if (t == 2) {
      dt = svdt;
    }
    if(non_block_comm)
    {
        //printf("---------------------------- Waiting for The Gathering ----------------------------\n");
        MPI_Waitall(48-12*((my_rank<NbCol)+(my_rank>=NbCol*(NbLi-1))+(my_rank%NbCol==0)+((my_rank+1)%NbCol==0)),reqs,stats) ;
    }
  }
  if(my_rank==0)
  if (file_export) {
  	printf("finalize_export\n");
    finalize_export(file);
  }

  free(hphy_send);
  free(hphy_recv);

  free(vphy_recv);
  free(vphy_send);

  free(uphy_send);
  free(uphy_recv);

  free(hfil_recv);
  free(hfil_send);

  free(ufil_recv);
  free(ufil_send);

  free(vfil_recv);
  free(vfil_send);

  free(hfil_buff_send);
  free(hfil_buff_recv);
}


void forward(void) {
  FILE *file = NULL;
  double svdt = 0.;
  int t = 0,k=0;
  MPI_Status status;
  MPI_Status stats[24];
  MPI_Request reqs[24];

  if(my_rank==0)
  {
	  if (file_export) {
	    //printf("1\n");
	    file = create_file();
	    //printf("2\n");
	    export_step(file, t);
	    //printf("3\n");
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

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < size_x/NP; i++) {
        for (int j = 0; j < size_y; j++) {
            if(my_rank==0)
            {
              hPhy_forward_vect(t,i,j);  
              uPhy_forward_vect(t,i,j);
              vPhy_forward_vect(t,i,j);
              hFil_forward_vect(t,i,j);
              uFil_forward_vect(t,i,j);
              vFil_forward_vect(t,i,j);
            }
            else
            {

              hPhy_forward_vect(t,i+1,j);  
              uPhy_forward_vect(t,i+1,j);
              vPhy_forward_vect(t,i+1,j);
              hFil_forward_vect(t,i+1,j);
              uFil_forward_vect(t,i+1,j);
              vFil_forward_vect(t,i+1,j);

            }
        }
    }

    //for(k=0;k<1;k++)
    {
        if(non_block_comm)
        {
            //printf("t=%d,k=%d:Communication non bloquante\n",t,k);
            if(my_rank!=0)
            {
                 MPI_Isend(&HPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_LAST_H_P,MPI_COMM_WORLD,&reqs[0]);
                 MPI_Irecv(&HPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_P, MPI_COMM_WORLD,&reqs[1]);

                 MPI_Isend(&UPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_LAST_U_P,MPI_COMM_WORLD,&reqs[2]);
                 MPI_Irecv(&UPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_U_P, MPI_COMM_WORLD,&reqs[3]);

                 MPI_Isend(&VPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_LAST_V_P,MPI_COMM_WORLD,&reqs[4]);
                 MPI_Irecv(&VPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_V_P, MPI_COMM_WORLD,&reqs[5]);

                 MPI_Isend(&HFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_LAST_H_F,MPI_COMM_WORLD,&reqs[6]);
                 MPI_Irecv(&HFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_F, MPI_COMM_WORLD,&reqs[7]);

                 MPI_Isend(&UFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_LAST_U_F,MPI_COMM_WORLD,&reqs[8]);
                 MPI_Irecv(&UFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_U_F, MPI_COMM_WORLD,&reqs[9]);

                 MPI_Isend(&VFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_LAST_V_F,MPI_COMM_WORLD,&reqs[10]);
                 MPI_Irecv(&VFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_V_F, MPI_COMM_WORLD,&reqs[11]);



            }
            if(my_rank!=NP-1)
            {
                 MPI_Isend(&HPHY_LOCAL(t + k,local_size_x-2,0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST_H_P
                           ,MPI_COMM_WORLD,&reqs[12-12*(my_rank==0)]);
                 MPI_Irecv(&HPHY_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_H_P
                           ,MPI_COMM_WORLD,&reqs[13-12*(my_rank==0)]);

                 MPI_Isend(&UPHY_LOCAL(t + k,local_size_x-2,0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST_U_P
                           ,MPI_COMM_WORLD,&reqs[14-12*(my_rank==0)]);
                 MPI_Irecv(&UPHY_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_U_P
                           ,MPI_COMM_WORLD,&reqs[15-12*(my_rank==0)]);

                 MPI_Isend(&VPHY_LOCAL(t + k,local_size_x-2,0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST_V_P
                           ,MPI_COMM_WORLD,&reqs[16-12*(my_rank==0)]);
                 MPI_Irecv(&VPHY_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_V_P
                           ,MPI_COMM_WORLD,&reqs[17-12*(my_rank==0)]);

                 MPI_Isend(&HFIL_LOCAL(t + k,local_size_x-2,0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST_H_F
                           ,MPI_COMM_WORLD,&reqs[18-12*(my_rank==0)]);
                 MPI_Irecv(&HFIL_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_H_F
                           ,MPI_COMM_WORLD,&reqs[19-12*(my_rank==0)]);

                 MPI_Isend(&UFIL_LOCAL(t + k,local_size_x-2,0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST_U_F
                           ,MPI_COMM_WORLD,&reqs[20-12*(my_rank==0)]);
                 MPI_Irecv(&UFIL_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_U_F
                           ,MPI_COMM_WORLD,&reqs[21-12*(my_rank==0)]);

                 MPI_Isend(&VFIL_LOCAL(t + k,local_size_x-2,0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST_V_F
                           ,MPI_COMM_WORLD,&reqs[22-12*(my_rank==0)]);
                 MPI_Irecv(&VFIL_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_V_F
                           ,MPI_COMM_WORLD,&reqs[23-12*(my_rank==0)]);


            }
        }
        else
        {
            if(my_rank!=0)
            {
                 MPI_Sendrecv(&HPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_H_P
                ,&HPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&UPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_U_P
                ,&UPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_U_P, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&VPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_V_P
                ,&VPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_V_P, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&UFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_U_F
                ,&UFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_U_F, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&VFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_V_F
                ,&VFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_V_F, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&HFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_H_F
                ,&HFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_F, MPI_COMM_WORLD,&status);

            }
            if(my_rank!=NP-1)
            {
                 MPI_Sendrecv(&HPHY_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_H_P
                ,&HPHY_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_H_P, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&UPHY_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_U_P
                ,&UPHY_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_U_P, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&VPHY_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_V_P
                ,&VPHY_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_V_P, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&UFIL_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_U_F
                ,&UFIL_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_U_F, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&VFIL_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_V_F
                ,&VFIL_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_V_F, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&HFIL_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_H_F
                ,&HFIL_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_H_F, MPI_COMM_WORLD,&status);

            }
        }
    }

    if(NP>1)
    {
        //printf("---------------------------- Magic The Gathering ----------------------------\n");
        MPI_Gather(&HFIL_LOCAL(t,(my_rank!=0), 0),size_y*size_x/NP,MPI_DOUBLE,
                   &HFIL(t, 0, 0),size_y*size_x/NP,MPI_DOUBLE,0,MPI_COMM_WORLD);
        //printf("---------------------------- End of Gathering  ----------------------------\n");
    }
	if(my_rank==0)
	{
	    if (file_export) {
	      export_step(file, t);
	    }
	    //printf("export_step\n");
	}
    if (t == 2) {
      dt = svdt;
    }
    if(non_block_comm)
    {
        //printf("---------------------------- Waiting for The Gathering ----------------------------\n");
        MPI_Waitall(24-12*(my_rank==0||my_rank==NP-1),reqs,stats) ;
    }
  }
  if(my_rank==0)
  if (file_export) {
  	//printf("finalize_export\n");
    finalize_export(file);
  }
}



void forward_parallel_io(void)
 {
  MPI_File file;// = NULL;
  double svdt = 0.;
  int t = 0,k=0;
  MPI_Status status;
  MPI_Status stats[24];
  MPI_Request reqs[24] ;

  //printf("---------MPI_IO---------\n");
  //if(my_rank==0)
  {
	  if (file_export) {
	    //printf("P#%d-create file\n",my_rank);
	    create_file_mpi(&file);
	    //printf("P#%d-export file t=%d\n",my_rank,t);
	    if(non_block_pararel_IO)
        {
            export_step_mpi_begin(&file,t);
        }
        else
            export_step_mpi(&file, t);

	    //printf("P#%d-file exported t=%d\n",my_rank,t);
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
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < size_x/NP; i++) {
        for (int j = 0; j < size_y; j++) {
          if(my_rank==0)
          {
            hPhy_forward_vect(t,i,j);  
            uPhy_forward_vect(t,i,j);
            vPhy_forward_vect(t,i,j);
            hFil_forward_vect(t,i,j);
            uFil_forward_vect(t,i,j);
            vFil_forward_vect(t,i,j);
          }
          else
          {

            hPhy_forward_vect(t,i+1,j);  
            uPhy_forward_vect(t,i+1,j);
            vPhy_forward_vect(t,i+1,j);
            hFil_forward_vect(t,i+1,j);
            uFil_forward_vect(t,i+1,j);
            vFil_forward_vect(t,i+1,j);

          }
      }
    }


    //for(k=0;k<1;k++)
    {
        if(non_block_comm)
        {
            //printf("t=%d,k=%d:Communication non bloquante\n",t,k);
            if(my_rank!=0)
            {
                 MPI_Isend(&HPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_LAST_H_P,MPI_COMM_WORLD,&reqs[0]);
                 MPI_Irecv(&HPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_P, MPI_COMM_WORLD,&reqs[1]);

                 MPI_Isend(&UPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_LAST_U_P,MPI_COMM_WORLD,&reqs[2]);
                 MPI_Irecv(&UPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_U_P, MPI_COMM_WORLD,&reqs[3]);

                 MPI_Isend(&VPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_LAST_V_P,MPI_COMM_WORLD,&reqs[4]);
                 MPI_Irecv(&VPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_V_P, MPI_COMM_WORLD,&reqs[5]);

                 MPI_Isend(&HFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_LAST_H_F,MPI_COMM_WORLD,&reqs[6]);
                 MPI_Irecv(&HFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_F, MPI_COMM_WORLD,&reqs[7]);

                 MPI_Isend(&UFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_LAST_U_F,MPI_COMM_WORLD,&reqs[8]);
                 MPI_Irecv(&UFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_U_F, MPI_COMM_WORLD,&reqs[9]);

                 MPI_Isend(&VFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_LAST_V_F,MPI_COMM_WORLD,&reqs[10]);
                 MPI_Irecv(&VFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_V_F, MPI_COMM_WORLD,&reqs[11]);



            }
            if(my_rank!=NP-1)
            {
                 MPI_Isend(&HPHY_LOCAL(t + k,local_size_x-2,0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST_H_P
                           ,MPI_COMM_WORLD,&reqs[12-12*(my_rank==0)]);
                 MPI_Irecv(&HPHY_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_H_P
                           ,MPI_COMM_WORLD,&reqs[13-12*(my_rank==0)]);

                 MPI_Isend(&UPHY_LOCAL(t + k,local_size_x-2,0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST_U_P
                           ,MPI_COMM_WORLD,&reqs[14-12*(my_rank==0)]);
                 MPI_Irecv(&UPHY_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_U_P
                           ,MPI_COMM_WORLD,&reqs[15-12*(my_rank==0)]);

                 MPI_Isend(&VPHY_LOCAL(t + k,local_size_x-2,0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST_V_P
                           ,MPI_COMM_WORLD,&reqs[16-12*(my_rank==0)]);
                 MPI_Irecv(&VPHY_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_V_P
                           ,MPI_COMM_WORLD,&reqs[17-12*(my_rank==0)]);

                 MPI_Isend(&HFIL_LOCAL(t + k,local_size_x-2,0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST_H_F
                           ,MPI_COMM_WORLD,&reqs[18-12*(my_rank==0)]);
                 MPI_Irecv(&HFIL_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_H_F
                           ,MPI_COMM_WORLD,&reqs[19-12*(my_rank==0)]);

                 MPI_Isend(&UFIL_LOCAL(t + k,local_size_x-2,0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST_U_F
                           ,MPI_COMM_WORLD,&reqs[20-12*(my_rank==0)]);
                 MPI_Irecv(&UFIL_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_U_F
                           ,MPI_COMM_WORLD,&reqs[21-12*(my_rank==0)]);

                 MPI_Isend(&VFIL_LOCAL(t + k,local_size_x-2,0),size_y, MPI_DOUBLE,my_rank+1,TAG_FIRST_V_F
                           ,MPI_COMM_WORLD,&reqs[22-12*(my_rank==0)]);
                 MPI_Irecv(&VFIL_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_V_F
                           ,MPI_COMM_WORLD,&reqs[23-12*(my_rank==0)]);


            }
        }
        else
        {
            if(my_rank!=0)
            {
                 MPI_Sendrecv(&HPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_H_P
                ,&HPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&UPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_U_P
                ,&UPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_U_P, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&VPHY_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_V_P
                ,&VPHY_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_V_P, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&UFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_U_F
                ,&UFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_U_F, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&VFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_V_F
                ,&VFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_V_F, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&HFIL_LOCAL(t + k,1, 0),size_y, MPI_DOUBLE, my_rank-1,TAG_LAST_H_F
                ,&HFIL_LOCAL(t + k,0, 0),size_y, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_F, MPI_COMM_WORLD,&status);

            }
            if(my_rank!=NP-1)
            {
                 MPI_Sendrecv(&HPHY_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_H_P
                ,&HPHY_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_H_P, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&UPHY_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_U_P
                ,&UPHY_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_U_P, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&VPHY_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_V_P
                ,&VPHY_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_V_P, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&UFIL_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_U_F
                ,&UFIL_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_U_F, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&VFIL_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_V_F
                ,&VFIL_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_V_F, MPI_COMM_WORLD,&status);

                 MPI_Sendrecv(&HFIL_LOCAL(t + k,local_size_x-2, 0),size_y, MPI_DOUBLE, my_rank+1,TAG_FIRST_H_F
                ,&HFIL_LOCAL(t + k,local_size_x-1,0),size_y, MPI_DOUBLE,my_rank+1,TAG_LAST_H_F, MPI_COMM_WORLD,&status);

            }
        }
    }

	{
	    if (file_export) {
            //printf("P#%d-export file t=%d\n",my_rank,t);
            if(non_block_pararel_IO)
            {
                export_step_mpi_end(&file,t-1);
                export_step_mpi_begin(&file,t);
            }
            else
                export_step_mpi(&file, t);
            //printf("P#%d-file exported for t=%d\n",my_rank,t);
	    }
	    //printf("export_step\n");
	}
    if (t == 2) {
      dt = svdt;
    }
    if(non_block_comm)
    {
        //printf("---------------------------- Waiting for The Gathering ----------------------------\n");
        MPI_Waitall(24-12*(my_rank==0||my_rank==NP-1),reqs,stats) ;
    }
  }
  //if(my_rank==0)
  if (file_export) {
  	//printf("finalize_export\n");
  	export_step_mpi_end(&file,t-1);
    finalize_export_mpi(&file);
  }
}



