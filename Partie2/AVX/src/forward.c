#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>

double hFil_forward(int t, int i, int j) {
 //Phase d'initialisation du filtre
 //HPHY(t - 1, i, j) est encore nul
 if (t <= 2)
   return HPHY(t, i, j);
 return HPHY(t - 1, i, j) +
   alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

void hFil_forward_vect(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //HPHY(t - 1, i, j) est encore nul

  __m256d hPhy_vect2 = _mm256_load_pd(&HPHY(t, i, j*4));
  //__m256d res = _mm256_setzero_ps();
  // printf("Hfill_vect(%d,%d,%d)=%lf\n",t, i, j,HFIL(t, 128, 128));
  // printf("Hfill_norm(%d,%d,%d)=%lf\n",t, i, j,hFil_forward(t, 128, 128));
  // printf("res:before return:%lf\n",HFIL(t, i, j+1));
  // printf("res:before return:%lf\n",HFIL(t, i, j+2));
  // printf("res:before return:%lf\n",HFIL(t, i, j+3));
  // printf("alpha:%lf\n",alpha);


  if (t <= 2)
  {
  	_mm256_store_pd(&HFIL(t, i, j*4),hPhy_vect2);//return HPHY(t, i, j);
  	return;
  }

  __m256d hFil_vect = _mm256_load_pd(&HFIL(t - 1, i, j*4));
  __m256d hPhy_vect = _mm256_load_pd(&HPHY(t - 1, i, j*4));
  __m256d res = _mm256_mul_pd(_mm256_set1_pd(alpha),hFil_vect);//alpha * (HFIL(t - 1, i, j)

  //hPhy_vect = _mm256_mul_pd(_mm256_set1_pd(2),hPhy_vect);//2*HPHY(t - 1, i, j)
  res = _mm256_sub_pd(_mm256_add_pd(hPhy_vect,res),_mm256_mul_pd(_mm256_set1_pd(2.0),hPhy_vect));
  //HPHY(t - 1, i, j)+alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j)
  res = _mm256_add_pd(res,hPhy_vect2);//HPHY(t - 1, i, j)+alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
   _mm256_store_pd(&HFIL(t, i, j*4),res);
  return;
//  HPHY(t - 1, i, j) +
//    alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

double uFil_forward(int t, int i, int j) {
 //Phase d'initialisation du filtre
 //UPHY(t - 1, i, j) est encore nul

 if (t <= 2)
   return UPHY(t, i, j);
 return UPHY(t - 1, i, j) +
   alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

void uFil_forward_vect(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //UPHY(t - 1, i, j) est encore nul
  __m256d uPhy_vect = _mm256_load_pd(&UPHY(t, i, j*4));

  if (t <= 2)
  {
  	_mm256_store_pd(&UFIL(t, i, j*4),uPhy_vect);
  	return;
  }

  __m256d uPhy_vect2 = _mm256_load_pd(&UPHY(t-1, i, j*4));
  __m256d uFil_vect  = _mm256_load_pd(&UFIL(t-1, i, j*4));

  __m256d res = _mm256_mul_pd(_mm256_set1_pd(alpha),uFil_vect);//alpha * (UFIL(t - 1, i, j)
  res = _mm256_sub_pd(_mm256_add_pd(uPhy_vect2,res),_mm256_mul_pd(_mm256_set1_pd(2.0),uPhy_vect2));
  //UPHY(t - 1, i, j)+alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j)
  res = _mm256_add_pd(res,uPhy_vect);//HPHY(t - 1, i, j)+alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
  _mm256_store_pd(&UFIL(t, i, j*4),res);
  return;
}

// double vFil_forward(int t, int i, int j) {
//  //Phase d'initialisation du filtre
//  //VPHY(t - 1, i, j) est encore nul
//  if (t <= 2)
//    return VPHY(t, i, j);
//  return VPHY(t - 1, i, j) +
//    alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
// }

void vFil_forward_vect(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //VPHY(t - 1, i, j) est encore nul
  __m256d vPhy_vect = _mm256_load_pd(&VPHY(t, i, j*4));

  if (t <= 2)
  {
  	_mm256_store_pd(&VFIL(t, i, j*4),vPhy_vect);
  	return;//VPHY(t, i, j);
  }

  __m256d vPhy_vect2 = _mm256_load_pd(&VPHY(t-1, i, j*4));
  __m256d vFil_vect  = _mm256_load_pd(&VFIL(t-1, i, j*4));

  __m256d res = _mm256_mul_pd(_mm256_set1_pd(alpha),vFil_vect);//alpha * (VFIL(t - 1, i, j)
  res = _mm256_sub_pd(_mm256_add_pd(vPhy_vect2,res),_mm256_mul_pd(_mm256_set1_pd(2.0),vPhy_vect2));
  //VPHY(t - 1, i, j)+alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j)
  res = _mm256_add_pd(res,vPhy_vect);//VPHY(t - 1, i, j)+alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
  _mm256_store_pd(&VFIL(t, i, j*4),res);
  return;
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

//inline
void hPhy_forward_vect(int t, int i, int j) {
  __m256d c, d;

  if (i > 0)
    c = _mm256_load_pd(&UPHY(t - 1, i - 1, j*4));
  else
    c=_mm256_setzero_pd();//_mm256_set1_pd(0.0);

  if (j < size_y/4 - 1)
    d = _mm256_loadu_pd(&VPHY(t - 1, i, j*4 + 1));
  else
    d = _mm256_set_pd(VPHY(t - 1, i, j*4 + 1),VPHY(t - 1, i, j*4 + 2),VPHY(t - 1, i, j*4 + 3),0.0);

  __m256d hFil_vect = _mm256_load_pd(&HFIL(t-1, i, j*4));
  __m256d uPhy_vect = _mm256_load_pd(&UPHY(t-1, i, j*4));
  __m256d vPhy_vect = _mm256_load_pd(&VPHY(t-1, i, j*4));

  __m256d res=_mm256_sub_pd(uPhy_vect,c);//(UPHY(t - 1, i, j) - c
  __m256d res2=_mm256_sub_pd(d,vPhy_vect);//d - VPHY(t - 1, i, j)
  //printf("dy:%lf\n",1.0/ dy);
  res2 = _mm256_mul_pd(_mm256_set1_pd(1.0/dy),res2);//(d - VPHY(t - 1, i, j)) / dy)
  res = _mm256_mul_pd(_mm256_set1_pd(dt*hmoy/ dx),res);
  //dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx
  //printf("dt*hmoy/ dx:%lf\n",dt*hmoy/ dx);
  res = _mm256_add_pd(_mm256_sub_pd(hFil_vect,res),res2);//HFIL(t - 1, i, j) -dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx +(d - VPHY(t - 1, i, j)) / dy)
  _mm256_store_pd(&HPHY(t, i, j*4),res);
  return;
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

void uPhy_forward_vect(int t, int i, int j) {
  __m256d b, e, f, g;

  if (i == size_x - 1)
  {
  	_mm256_store_pd(&UPHY(t, i, j*4),_mm256_setzero_pd());
  	return;
  }

  if (i < size_x - 1)
    b = _mm256_load_pd(&HPHY(t - 1, i + 1, j*4));
  else
    b=_mm256_setzero_pd();

  if (j < size_y/4 - 1)
    e = _mm256_loadu_pd(&VPHY(t - 1, i, j*4 + 1));
  else
    e=_mm256_set_pd(VPHY(t - 1, i, j*4 + 1),VPHY(t - 1, i, j*4 + 2),VPHY(t - 1, i, j*4 + 3),0.0);

  if (i < size_x - 1)
    f = _mm256_load_pd(&VPHY(t - 1, i + 1, j*4));
  else
    f = _mm256_setzero_pd();

  if (i < size_x - 1 && j < size_y/4 - 1)
    g = _mm256_loadu_pd(&VPHY(t - 1, i + 1, j*4 + 1));
  else
  {
    if(i < size_x -1)
        g = _mm256_set_pd(VPHY(t - 1, i + 1, j*4 + 1),
        VPHY(t - 1, i + 1, j*4 + 2),
        VPHY(t - 1, i + 1, j*4 + 3),0.0);
    else
        g = _mm256_setzero_pd();
  }
    //_mm256_set_pd(VPHY(t - 1, i + 1, j*4 + 1),VPHY(t - 1, i + 1, j*4 + 2),VPHY(t - 1, i + 1, j*4 + 3),0.0);

  __m256d uFil_vect=_mm256_load_pd(&UFIL(t - 1, i, j*4));
  __m256d hPhy_vect=_mm256_load_pd(&HPHY(t - 1, i, j*4));
  __m256d vPhy_vect=_mm256_load_pd(&VPHY(t - 1, i, j*4));

  __m256d res=_mm256_sub_pd(b,hPhy_vect);//b - HPHY(t - 1, i, j)
  __m256d res2=_mm256_add_pd(e,_mm256_add_pd(f,g));//e + f + g
  res2=_mm256_mul_pd(_mm256_set1_pd(pcor / 4. ),_mm256_add_pd(vPhy_vect,res2));
  //(pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g)
  res=_mm256_fmadd_pd(_mm256_set1_pd(-grav / dx),res,res2);
  //(-grav / dx) * (b - HPHY(t - 1, i, j)) + (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g)
  res=_mm256_sub_pd(res,_mm256_mul_pd(_mm256_set1_pd(dissip),uFil_vect));
  //(-grav / dx) * (b - HPHY(t - 1, i, j)) + (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) -
  //(dissip * UFIL(t - 1, i, j))
  res=_mm256_mul_pd(_mm256_set1_pd(dt),res);//dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) +
  //(pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) -
  //(dissip * UFIL(t - 1, i, j)))
  res=_mm256_add_pd(res,uFil_vect);
//  UFIL(t - 1, i, j) +
//    dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) +
//	  (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) -
//	  (dissip * UFIL(t - 1, i, j)))
  _mm256_store_pd(&UPHY(t, i, j*4),res);

  return;
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

void vPhy_forward_vect(int t, int i, int j) {
  __m256d c, d, e, f;

  if (j == 0)
  {
  	// _mm256_store_pd(&VPHY(t, i, j*4),_mm256_set_pd(0.0,vPhy_forward(t,i,1),vPhy_forward(t,i,2),vPhy_forward(t,i,3)));
  	_mm256_store_pd(&VPHY(t, i, j*4),_mm256_setzero_pd());
  	vPhy_forward_vect(t,i,j+1);
  	return;
  }

  if (j > 0)
    c = _mm256_loadu_pd(&HPHY(t - 1, i, j*4 - 1));
  else
    c = _mm256_set_pd(0.0,HPHY(t - 1, i, 0),HPHY(t - 1, i , 1),HPHY(t - 1, i, 2));

  if (i > 0 && j > 0)
    d = _mm256_loadu_pd(&UPHY(t - 1, i - 1, j*4 -1));
  else
  {
    if(i>0)
        d = _mm256_set_pd(0.0,UPHY(t - 1, i - 1, 0),UPHY(t - 1, i - 1 , 1),UPHY(t - 1, i - 1, 2));
    else
        d = _mm256_setzero_pd();
  }

  if (i > 0)
    e = _mm256_load_pd(&UPHY(t - 1, i - 1, j*4));
  else
    e = _mm256_setzero_pd();

  if (j > 0)
    f = _mm256_loadu_pd(&UPHY(t - 1, i, j*4 - 1));
  else
    f = _mm256_set_pd(0.0,UPHY(t - 1, i, 0),UPHY(t - 1, i , 1),UPHY(t - 1, i, 2));

  __m256d vFil_vect=_mm256_load_pd(&VFIL(t - 1, i, j*4));
  __m256d hPhy_vect=_mm256_load_pd(&HPHY(t - 1, i, j*4));
  __m256d uPhy_vect=_mm256_load_pd(&UPHY(t - 1, i, j*4));

  __m256d res=_mm256_sub_pd(hPhy_vect,c);//HPHY(t - 1, i, j) - c
  __m256d res2=_mm256_add_pd(d,_mm256_add_pd(e,f));//d + e + f
  res2=_mm256_mul_pd(_mm256_set1_pd(pcor / 4. ),_mm256_add_pd(uPhy_vect,res2));
  //(pcor / 4.) * (d + e + f + UPHY(t - 1, i, j))
  res=_mm256_fmsub_pd(_mm256_set1_pd(-grav / dy),res,res2);
  //(-grav / dy) * (HPHY(t - 1, i, j) - c) - (pcor / 4.) * (d + e + f + UPHY(t - 1, i, j))
  res=_mm256_sub_pd(res,_mm256_mul_pd(_mm256_set1_pd(dissip),vFil_vect));
  //(-grav / dy) * (HPHY(t - 1, i, j) - c) - (pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) -
  //(dissip * UFIL(t - 1, i, j))
  res=_mm256_mul_pd(_mm256_set1_pd(dt),res);
  // dt * ((-grav / dy) * (HPHY(t - 1, i, j) - c) -
  // (pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) -
  // (dissip * VFIL(t - 1, i, j)))
  res=_mm256_add_pd(res,vFil_vect);
  // VFIL(t - 1, i, j) +
  // dt * ((-grav / dy) * (HPHY(t - 1, i, j) - c) -
  // (pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) -
  // (dissip * VFIL(t - 1, i, j)));
  _mm256_store_pd(&VPHY(t, i, j*4),res);
  return;
}

void forward(void) {
  FILE *file = NULL;
  double svdt = 0.;
  int t = 0;
  //__m256d hFil_vect, uFil_vect, vFil_vect, hPhy_vect, uPhy_vect, vPhy_vect;

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
		// HPHY(t, i, j) = hPhy_forward(t, i, j);
  //       UPHY(t, i, j) = uPhy_forward(t, i, j);
  //       VPHY(t, i, j) = vPhy_forward(t, i, j);
  //       HFIL(t, i, j) = hFil_forward(t, i, j);
  //       UFIL(t, i, j) = uFil_forward(t, i, j);
  //       VFIL(t, i, j) = vFil_forward(t, i, j);

//	vFil_vect = _mm256_load_pd(&VFIL(t, i, j*4));
//	uFil_vect = _mm256_load_pd(&UFIL(t, i, j*4));
//	hFil_vect = _mm256_load_pd(&HFIL(t, i, j*4));
//	vPhy_vect = _mm256_load_pd(&VPHY(t, i, j*4));
//	uPhy_vect = _mm256_load_pd(&UPHY(t, i, j*4));
//	hPhy_vect = _mm256_load_pd(&HPHY(t, i, j*4));
//	//_mm256_store_ps(&VFIL(t, i, j*4),v3);
    //printf("HFIl:after return:%lf\n",hFil_vect[0]);
    //   	if(i<4 && j%4==0)
    //   	{
    //   		printf("----------------------------------------------\n");
    //   		printf("before uphy(%d,%d,%d)=%lf\n",t, i, j,HPHY(t, i,j));
    //   	}
    
    // if(i<4 && j%4==0 && j<16)
    // 	printf("expected after uphy(%d,%d,%d)=%lf\n",t, i, j,	hPhy_forward(t, i,j));
    hPhy_forward_vect(t,i,j);  
    uPhy_forward_vect(t,i,j);
    vPhy_forward_vect(t,i,j);
    hFil_forward_vect(t,i,j);
    uFil_forward_vect(t,i,j);
    vFil_forward_vect(t,i,j);
   //  if(i<4 && j<8)
  	// {
  	// 	printf("after hphy(%d,%d,%d)=%lf\n",t, i, j,HPHY(t, i,j));
  	// 	printf("after uphy(%d,%d,%d)=%lf\n",t, i, j,UPHY(t, i,j));
  	// 	printf("after vphy(%d,%d,%d)=%lf\n",t, i, j,VPHY(t, i,j));
  	// 	printf("after hFil(%d,%d,%d)=%lf\n",t, i, j,HFIL(t, i,j));
  	// 	printf("after ufil(%d,%d,%d)=%lf\n",t, i, j,UFIL(t, i,j));
  	// 	printf("after vfil(%d,%d,%d)=%lf\n",t, i, j,VFIL(t, i,j));
  	// 	printf("----------------------------------------------\n");

  	// }
    // _mm256_store_pd(&HFIL(t, i, j*4),hFil_vect);
    // _mm256_store_pd(&UFIL(t, i, j*4),uFil_vect);
    // _mm256_store_pd(&VFIL(t, i, j*4),vFil_vect);
    // _mm256_store_pd(&HPHY(t, i, j*4),hPhy_vect);
    // _mm256_store_pd(&UPHY(t, i, j*4),uPhy_vect);
    // _mm256_store_pd(&VPHY(t, i, j*4),vPhy_vect);

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
