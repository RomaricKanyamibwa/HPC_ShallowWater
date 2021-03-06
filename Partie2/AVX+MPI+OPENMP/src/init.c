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

  #pragma omp parallel private(tmp)
  for (int i = 0; i < size_x;  i++) {
    for (int j = 0; j < size_y; j++) {
    tmp=0;
	if(i<size_x/NP+1)
	{
        tmp=(size_x/NP*my_rank)/*-1*(i>0)*(my_rank!=0)*/;
        //printf("P#%d:tmp=%lf\n",my_rank,tmp);
	}else
	{
        break;
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
  gmx = size_x * dx / 2 ;
  gmy = size_y * dy / 2 ;
  gsx = 25000 ;
  gsy = 25000 ;

  #pragma omp parallel
  for (int i = 0; i < size_x/NbLi;  i++) {
    for (int j = 0; j < size_y/NbCol; j++) {

	HFIL_LOCAL(0, i+(my_rank>=NbCol), j+(my_rank%NbCol!=0)) = height *
	(exp(- pow(((i+(my_rank/NbCol)*size_x/NbLi) * dx - gmx) / gsx, 2) / 2.)) *
	(exp(- pow(((j+(my_rank%NbCol)*size_y/NbCol) * dy - gmy) / gsy, 2) / 2.)) ;
    }
  }
}
