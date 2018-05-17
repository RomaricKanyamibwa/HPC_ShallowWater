#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>

//double hFil_forward(int t, int i, int j) {
//  //Phase d'initialisation du filtre
//  //HPHY(t - 1, i, j) est encore nul
//  if (t <= 2)
//    return HPHY(t, i, j);
//  return HPHY(t - 1, i, j) +
//    alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
//}

__m256d hFil_forward_vect(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //HPHY(t - 1, i, j) est encore nul

  __m256d hPhy_vect2 = _mm256_load_pd(&HPHY(t, i, j*4));
  //__m256d res = _mm256_setzero_ps();


  if (t <= 2)
    return hPhy_vect2;

  __m256d hFil_vect = _mm256_load_pd(&HFIL(t - 1, i, j*4));
  __m256d hPhy_vect = _mm256_load_pd(&HPHY(t - 1, i, j*4));
  __m256d res = _mm256_mul_pd(_mm256_set1_pd(alpha),hFil_vect);
  //alpha * (HFIL(t - 1, i, j)
  //hPhy_vect = _mm256_mul_pd(_mm256_set1_pd(2),hPhy_vect);//2*HPHY(t - 1, i, j)
  res = _mm256_sub_pd(res,hPhy_vect);//HPHY(t - 1, i, j)+alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j)
  res = _mm256_add_pd(res,hPhy_vect2);//HPHY(t - 1, i, j)+alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
  return res;
//  HPHY(t - 1, i, j) +
//    alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

//double uFil_forward(int t, int i, int j) {
//  //Phase d'initialisation du filtre
//  //UPHY(t - 1, i, j) est encore nul
//
//  if (t <= 2)
//    return UPHY(t, i, j);
//  return UPHY(t - 1, i, j) +
//    alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
//}

__m256d uFil_forward_vect(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //UPHY(t - 1, i, j) est encore nul
  __m256d uPhy_vect = _mm256_load_pd(&UPHY(t, i, j*4));

  if (t <= 2)
    return uPhy_vect;

  __m256d uPhy_vect2 = _mm256_load_pd(&UPHY(t-1, i, j*4));
  __m256d uFil_vect  = _mm256_load_pd(&UFIL(t-1, i, j*4));

  __m256d res = _mm256_mul_pd(_mm256_set1_pd(alpha),uFil_vect);//alpha * (UFIL(t - 1, i, j)
  res = _mm256_sub_pd(res,uPhy_vect2);//UPHY(t - 1, i, j)+alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j)
  res = _mm256_add_pd(res,uPhy_vect);//HPHY(t - 1, i, j)+alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
  return res;
}

//double vFil_forward(int t, int i, int j) {
//  //Phase d'initialisation du filtre
//  //VPHY(t - 1, i, j) est encore nul
//  if (t <= 2)
//    return VPHY(t, i, j);
//  return VPHY(t - 1, i, j) +
//    alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
//}

__m256d vFil_forward_vect(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //VPHY(t - 1, i, j) est encore nul
  __m256d vPhy_vect = _mm256_load_pd(&VPHY(t, i, j*4));

  if (t <= 2)
    return vPhy_vect;//VPHY(t, i, j);

  __m256d vPhy_vect2 = _mm256_load_pd(&VPHY(t-1, i, j*4));
  __m256d vFil_vect  = _mm256_load_pd(&VFIL(t-1, i, j*4));

  __m256d res = _mm256_mul_pd(_mm256_set1_pd(alpha),vFil_vect);//alpha * (VFIL(t - 1, i, j)
  res = _mm256_sub_pd(res,vPhy_vect2);//VPHY(t - 1, i, j)+alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j)
  res = _mm256_add_pd(res,vPhy_vect);//VPHY(t - 1, i, j)+alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
  return res;
}

//double hPhy_forward(int t, int i, int j) {
//  double c, d;
//
//  c = 0.;
//  if (i > 0)
//    c = UPHY(t - 1, i - 1, j);
//
//  d = 0.;
//  if (j < size_y - 1)
//    d = VPHY(t - 1, i, j + 1);
//
//  return HFIL(t - 1, i, j) -
//    dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx +
//		 (d - VPHY(t - 1, i, j)) / dy);
//}

//inline
__m256d hPhy_forward_vect(int t, int i, int j) {
  __m256d c, d;

  if (i > 0)
    c = _mm256_load_pd(&UPHY(t - 1, i - 1, j*4));
  else
    c=_mm256_set1_pd(0.0);

  if (j < size_y/4)
    d = _mm256_loadu_pd(&VPHY(t - 1, i, j*4 + 1));
  else
    d = _mm256_set_pd(VPHY(t - 1, i, j*4 + 1),VPHY(t - 1, i, j*4 + 2),VPHY(t - 1, i, j*4 + 3),0.0);

  __m256d hFil_vect = _mm256_load_pd(&UFIL(t-1, i, j*4));
  __m256d uPhy_vect = _mm256_load_pd(&UPHY(t-1, i, j*4));
  __m256d vPhy_vect = _mm256_load_pd(&VPHY(t-1, i, j*4));

  __m256d res=_mm256_sub_pd(uPhy_vect,c);//(UPHY(t - 1, i, j) - c
  __m256d res2=_mm256_sub_pd(d,vPhy_vect);//d - VPHY(t - 1, i, j)
  //printf("dy:%lf\n",1.0/ dy);
  res2 = _mm256_mul_pd(_mm256_set1_pd(1.0/ dy),res2);//(d - VPHY(t - 1, i, j)) / dy)
  res = _mm256_fmadd_pd(_mm256_set1_pd(dt*hmoy/ dx),res,res2);//dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx +(d - VPHY(t - 1, i, j)) / dy);
  //printf("dt*hmoy/ dx:%lf\n",dt*hmoy/ dx);
  res = _mm256_add_pd(hFil_vect,res);//HFIL(t - 1, i, j) -dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx +(d - VPHY(t - 1, i, j)) / dy)
  return res;
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

__mm256d uPhy_forward(int t, int i, int j) {
  __mm256d b, e, f, g;

  if (i == size_x - 1)
    return _mm256_setzero_ps();

  if (i < size_x - 1)
    b = HPHY(t - 1, i + 1, j);
  else
    b=_mm256_setzero_ps();

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
  int t = 0;
  __m256d hFil_vect, uFil_vect, vFil_vect, hPhy_vect, uPhy_vect, vPhy_vect;

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

    for (int i = 0; i < size_x; i++) {
      for (int j = 0; j < size_y/4; j++) {
	HPHY(t, i, j) = hPhy_forward(t, i, j);
	UPHY(t, i, j) = uPhy_forward(t, i, j);
	VPHY(t, i, j) = vPhy_forward(t, i, j);
//	HFIL(t, i, j) = hFil_forward(t, i, j);
//	UFIL(t, i, j) = uFil_forward(t, i, j);
//	VFIL(t, i, j) = vFil_forward(t, i, j);

//	vFil_vect = _mm256_load_pd(&VFIL(t, i, j*4));
//	uFil_vect = _mm256_load_pd(&UFIL(t, i, j*4));
//	hFil_vect = _mm256_load_pd(&HFIL(t, i, j*4));
//	vPhy_vect = _mm256_load_pd(&VPHY(t, i, j*4));
//	uPhy_vect = _mm256_load_pd(&UPHY(t, i, j*4));
//	hPhy_vect = _mm256_load_pd(&HPHY(t, i, j*4));
//	//_mm256_store_ps(&VFIL(t, i, j*4),v3);
	__m256d res=_mm256_mul_pd(vFil_vect,vPhy_vect);
	res=_mm256_mul_pd(uFil_vect,uPhy_vect);
	//res=_mm256_mul_pd(hFil_vect,hPhy_vect);
	res=_mm256_mul_pd(res,res);
    //hPhy_vect=hPhy_forward_vect(t,i,j);
    hFil_vect=hFil_forward_vect(t,i,j);
    uFil_vect=uFil_forward_vect(t,i,j);
    vFil_vect=vFil_forward_vect(t,i,j);
    hPhy_vect=hPhy_forward_vect(t,i,j);
    _mm256_store_pd(&HFIL(t, i, j*4),hFil_vect);
    _mm256_store_pd(&UFIL(t, i, j*4),uFil_vect);
    _mm256_store_pd(&VFIL(t, i, j*4),vFil_vect);
    _mm256_store_pd(&HPHY(t, i, j*4),hPhy_vect);

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
