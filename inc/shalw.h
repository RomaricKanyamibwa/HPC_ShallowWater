#include <string>
#include <mpi.h>

extern double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
extern int size_x/**Numero de lignes**/, size_y, nb_steps,local_size_x, local_size_y/*bloc*/;
extern int NbCol, NbLi, ligne_colonne /*bloc*/; //nombre de lignes et colonnes pour le par bloc
extern double *hFil_local, *uFil_local, *vFil_local, *hPhy_local, *uPhy_local, *vPhy_local;
extern double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
extern bool file_export;
extern int my_rank,NP/*Nombre de processeur*/,root,NP_temp/*bloc*/;
extern std::string export_path;

#define HFIL(t, i, j) hFil[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define UFIL(t, i, j) uFil[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define VFIL(t, i, j) vFil[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define HPHY(t, i, j) hPhy[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define UPHY(t, i, j) uPhy[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define VPHY(t, i, j) vPhy[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
//Macro Local
//ATTENTION OPTI LOCAL_SIZEX
#define HFIL_LOCAL(t, i, j) hFil_local[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define UFIL_LOCAL(t, i, j) uFil_local[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define VFIL_LOCAL(t, i, j) vFil_local[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define HPHY_LOCAL(t, i, j) hPhy_local[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define UPHY_LOCAL(t, i, j) uPhy_local[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define VPHY_LOCAL(t, i, j) vPhy_local[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
