#include <math.h>
#include <shalw.h>

void gauss_init_loc(int rank,int p);

void gauss_init(void) {
  double gmx, gmy, gsx, gsy;
  double tmp=0.0;

  gmx = size_x * dx / 2 ;
  gmy = size_y * dy / 2 ;
  gsx = 25000 ;
  gsy = 25000 ;

  for (int i = 0; i < size_x;  i++) {
    for (int j = 0; j < size_y; j++) {
      HFIL(0, i, j) = height *
	(exp(- pow((i * dx - gmx) / gsx, 2) / 2.)) *
	(exp(- pow((j * dy - gmy) / gsy, 2) / 2.)) ;
	//ATTENTION OPTI LOCAL_SIZEX
	if(i<size_x/NP+1)
	{
        tmp=(size_x/NP*my_rank)/*-1*(i>0)*(my_rank!=0)*/;
        //printf("P#%d:tmp=%lf\n",my_rank,tmp);
	}
	HFIL_LOCAL(0, i, j) = height *
	(exp(- pow(((i+tmp) * dx - gmx) / gsx, 2) / 2.)) *
	(exp(- pow((j * dy - gmy) / gsy, 2) / 2.)) ;
    //gauss_init_loc(my_rank,NP);
    }
  }
}

void gauss_init_bloc(void) {
  double gmx, gmy, gsx, gsy;
  int tmpx=(size_x/NbLi*my_rank/NbCol);
  int tmpy=(size_y/NbCol*my_rank%NbCol);

  gmx = size_x * dx / 2 ;
  gmy = size_y * dy / 2 ;
  gsx = 25000 ;
  gsy = 25000 ;

  for (int i = 0; i < size_x/NbLi;  i++) {
    for (int j = 0; j < size_y/NbCol; j++) {
//      HFIL(0, i, j) = height *
//	(exp(- pow((i * dx - gmx) / gsx, 2) / 2.)) *
//	(exp(- pow((j * dy - gmy) / gsy, 2) / 2.)) ;

	HFIL_LOCAL(0, i, j) = height *
	(exp(- pow(((i+100) * dx - gmx) / gsx, 2) / 2.)) *
	(exp(- pow(((j+100) * dy - gmy) / gsy, 2) / 2.)) ;
	//ATTENTION OPTI LOCAL_SIZEX
//	if(i<size_x/NbLi)
//	{
//        tmpx=(size_x/NbLi*my_rank/NbCol)/*-1*(i>0)*(my_rank!=0)*/;
        //printf("P#%d:tmp=%lf\n",my_rank,tmp);
//	}
//	else
//        tmpx=0.0;
//	if(j<size_y/NbCol)
//	{
//        tmpy=(size_y/NbCol*my_rank%NbCol)/*-1*(i>0)*(my_rank!=0)*/;
        //printf("P#%d:tmp=%lf\n",my_rank,tmp);
//	}
//	else
//        tmpy=0.0;
//	HFIL_LOCAL(0, i/*+(my_rank>=NbCol)*/, j/*+(my_rank%NbCol!=0)*/) = height *
//	(exp(- pow(((i+tmpy) * dx - gmx) / gsx, 2) / 2.)) *
//	(exp(- pow(((j+tmpx) * dy - gmy) / gsy, 2) / 2.)) ;
    }
  }
  tmpx++;
  tmpy++;
}

void gauss_init_loc(int rank,int p) {
  double gmx, gmy, gsx, gsy;

  gmx = size_x * dx / 2 ;
  gmy = size_y * dy / 2 ;
  gsx = 25000 ;
  gsy = 25000 ;

  int k = size_x/p; // nb element

  int z = rank*k;   //debut de  i

  //printf("rank:%d ----- z:%d -------- z+k:%d\n",rank,z,z+k);

  for (int i = z; i < z+k;  i++) {
    for (int j = 0; j < size_y; j++) {

    hFil_local[ (j) + (i) * size_y -rank*k*size_y +(rank!=0)*size_y] = height *
  (exp(- pow((i * dx - gmx) / gsx, 2) / 2.)) *
  (exp(- pow((j * dy - gmy) / gsy, 2) / 2.)) ;

  //printf("rank %d   val init :%f   %d\n",rank,hFil_loc[ (j) + (i) * size_y -rank*k*size_y +(rank!=0)*size_y], (j) + (i) * size_y);
    }
  }
}
