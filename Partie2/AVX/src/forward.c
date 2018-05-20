#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>
#include <immintrin.h>

void hFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //HPHY(t - 1, i, j) est encore nul
  __m256d a, b, c,r;
  if (t <= 2){
  	a = _mm256_load_pd(&HPHY(t, i, j*4));
  	_mm256_store_pd(&HFIL(t,i,j*4), a);

  }
  else{
  	a = _mm256_load_pd(&HPHY(t - 1, i, j*4));
	b = _mm256_load_pd(&HPHY(t, i, j*4));
	c = _mm256_load_pd(&HFIL(t - 1, i, j*4));
	r = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(alpha),c),_mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(-2),a),b));
    _mm256_store_pd(&HFIL(t,i,j*4), r);
  }
}


void uFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //UPHY(t - 1, i, j) est encore nul
  __m256d a, b, c,r;
  if(t<=2){
  		a = _mm256_load_pd(&UPHY(t, i, j*4));
  	   _mm256_store_pd(&UFIL(t,i,j*4), a);
  }else{

  		a = _mm256_load_pd(&UPHY(t - 1, i, j*4));
		b = _mm256_load_pd(&UPHY(t, i, j*4));
		c = _mm256_load_pd(&UFIL(t - 1, i, j*4));
		r = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(alpha),c),_mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(-2),a),b));
    	_mm256_store_pd(&UFIL(t,i,j*4), r);
  }

 // if (t <= 2)
 //  return UPHY(t, i, j);
 // return UPHY(t - 1, i, j) +
 //   alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

void vFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //VPHY(t - 1, i, j) est encore nul
  __m256d a, b, c,r;

  if(t<=2){
  		a = _mm256_load_pd(&VPHY(t, i, j*4));
  	   _mm256_store_pd(&VFIL(t,i,j*4), a);
  }else{

  		a = _mm256_load_pd(&VPHY(t - 1, i, j*4));
		b = _mm256_load_pd(&VPHY(t, i, j*4));
		c = _mm256_load_pd(&VFIL(t - 1, i, j*4));
		r = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(alpha),c),_mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(-2),a),b));;
    	_mm256_store_pd(&VFIL(t,i,j*4), r);
  }


  //if (t <= 2)
  //  return VPHY(t, i, j);
  //return VPHY(t - 1, i, j) +
  //  alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
}

void hPhy_forward(int t, int i, int j) {
  //double c, d;
  __m256d hf,up,vp, c, d, r,r1;

 // c = 0.;

 c =_mm256_set1_pd(0);

  if (i > 0)
    //c = UPHY(t - 1, i - 1, j);
    c = _mm256_load_pd(&UPHY(t - 1, i-1, j*4));
	

  //d = 0.;
d =  _mm256_set1_pd(0);
  if (j < size_y/4-1)
    //d = VPHY(t - 1, i, j + 1);
    d = _mm256_load_pd(&VPHY(t - 1, i, (j+1)*4));



hf =  _mm256_load_pd(&HFIL(t - 1, i, j*4));
up =  _mm256_load_pd(&UPHY(t - 1, i, j*4));
vp =  _mm256_load_pd(&VPHY(t - 1, i, j*4));

r = _mm256_add_pd(up, _mm256_div_pd( _mm256_mul_pd(_mm256_set1_pd(-1),c),_mm256_set1_pd(dx)));
r = _mm256_mul_pd(_mm256_set1_pd(hmoy),_mm256_mul_pd(r,_mm256_set1_pd(dt)));
r1 = _mm256_add_pd(d,_mm256_div_pd( _mm256_mul_pd(_mm256_set1_pd(-1),vp),_mm256_set1_pd(dy)));
r = _mm256_add_pd(r,hf);
r = _mm256_add_pd(r,r1);
_mm256_store_pd(&HPHY(t,i,j*4), r);

//  return HFIL(t - 1, i, j) -
//    dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx +
//		 (d - VPHY(t - 1, i, j)) / dy);
}

void uPhy_forward(int t, int i, int j) {
  //double b, e, f, g;
  __m256d uf,hp,vp,r,r1, b,e,f,g;
  if (i == size_x - 1){
  	b =_mm256_set1_pd(0);
  	_mm256_store_pd(&UPHY(t,i,j*4), b);

  }
    

  b =_mm256_set1_pd(0);
  if (i < size_x - 1)
    //b = HPHY(t - 1, i + 1, j);
    b = _mm256_load_pd(&HPHY(t - 1, i+1, j*4));

 	e =_mm256_set1_pd(0);
  if (j < size_y - 1)
    //e = VPHY(t - 1, i, j + 1);
	e = _mm256_load_pd(&VPHY(t - 1, i, (j+1)*4));

  f =_mm256_set1_pd(0);
  if (i < size_x - 1)
    //f = VPHY(t - 1, i + 1, j);
	f = _mm256_load_pd(&VPHY(t - 1, i+1, (j)*4));

  g =_mm256_set1_pd(0);
  if (i < size_x - 1 && j < size_y - 1)
    //g = VPHY(t - 1, i + 1, j + 1);
    g = _mm256_load_pd(&VPHY(t - 1, i+1, (j+1)*4));


  uf =  _mm256_load_pd(&HFIL(t - 1, i, j*4));
  hp =  _mm256_load_pd(&UPHY(t - 1, i, j*4));
  vp =  _mm256_load_pd(&VPHY(t - 1, i, j*4));

  r = _mm256_mul_pd(_mm256_set1_pd(dt),_mm256_div_pd(_mm256_set1_pd(-grav),_mm256_set1_pd(dx)));
  r = _mm256_mul_pd(r,_mm256_add_pd(b,hp));


  r1 =  _mm256_add_pd(e,f);
  r1 =  _mm256_add_pd(r1,g);
  r1 =   _mm256_add_pd(r1,vp);//(pcor / 4.) *
  r1 = _mm256_mul_pd(_mm256_div_pd(_mm256_set1_pd(pcor),_mm256_set1_pd(4.0)) ,_mm256_add_pd(r1,vp));
  r1 =   _mm256_add_pd(r1,_mm256_mul_pd(_mm256_set1_pd(dissip),uf));
  r =    _mm256_add_pd(r,_mm256_mul_pd(_mm256_set1_pd(-1),r1));
  r =    _mm256_add_pd(uf,r);

  _mm256_store_pd(&UPHY(t,i,j*4), r);



  //return UFIL(t - 1, i, j) +
  //  dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) +
	//  (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) -
	//  (dissip * UFIL(t - 1, i, j)));
}

void vPhy_forward(int t, int i, int j) {
  
   __m256d vf,hp,up,r,r1, c, d, e, f,r2;

  c =_mm256_set1_pd(0);
  if (j == 0)
  	 _mm256_store_pd(&VPHY(t,i,j*4), c);

  //c = 0.;
  c =_mm256_set1_pd(0);
  if (j > 0)
    //c = HPHY(t - 1, i, j - 1);
	c = _mm256_load_pd(&HPHY(t - 1, i-1, (j-1)*4));

  d =_mm256_set1_pd(0);
  if (i > 0 && j > 0)
    //d = UPHY(t - 1, i -1, j -1);
	d = _mm256_load_pd(&UPHY(t - 1, i-1, (j-1)*4));

   e =_mm256_set1_pd(0);
  if (i > 0)
    //e = UPHY(t - 1, i - 1, j);
	e = _mm256_load_pd(&UPHY(t - 1, i-1, j*4));

   f =_mm256_set1_pd(0);
  if (j > 0)
    //f = UPHY(t - 1, i, j - 1);
	f = _mm256_load_pd(&UPHY(t - 1, i, (j-1)*4));

  
  vf =  _mm256_load_pd(&HFIL(t - 1, i, j*4));
  hp =  _mm256_load_pd(&UPHY(t - 1, i, j*4));
  up =  _mm256_load_pd(&VPHY(t - 1, i, j*4));

  r =  _mm256_add_pd(d,e);
  r1 =  _mm256_add_pd(f,up);
  r =  _mm256_add_pd(r,r1);
  r =  _mm256_add_pd(r,-dissip *vf); // (pcor / 4.) *
  r2 = _mm256_div_pd(_mm256_set1_pd(pcor),_mm256_set1_pd(4.0));
  r = _mm256_mul_pd(r,r2);

  r1 = _mm256_mul_pd(_mm256_set1_pd(dt),_mm256_div_pd(_mm256_set1_pd(-grav),_mm256_set1_pd(dy)));
  r1 = _mm256_mul_pd(r1,_mm256_add_pd(hp,_mm256_mul_pd(_mm256_set1_pd(-1),c)));
  r =  _mm256_add_pd(r,r1);
  r =  _mm256_add_pd(r,vf);
  _mm256_store_pd(&VPHY(t,i,j*4), r);

  //return VFIL(t - 1, i, j) +
   // dt * ((-grav / dy) * (HPHY(t - 1, i, j) - c) -
	//  (pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) -
	//  (dissip * VFIL(t - 1, i, j)));
}

void forward(void) {
  FILE *file = NULL;
  double svdt = 0.;
  int t = 0;
  
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

     for (int i = 0; i < size_x; i++){ // inversion des boucles
       for (int j = 0; j < size_y/4; j++){
		hPhy_forward(t, i, j);
	 	uPhy_forward(t, i, j);
		vPhy_forward(t, i, j);
		hFil_forward(t, i, j);
		uFil_forward(t, i, j);
		vFil_forward(t, i, j);
		
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

