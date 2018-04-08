void gauss_init_bloc(void) {
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



void forward_bloc(void) {
  FILE *file = NULL;
  double svdt = 0.;
  int t = 0,k,i;
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
    	/* au dessus en dessous */
        if(my_rank>=NbCol) //on envoie celui du haut
        {
            mpi_ret_type = MPI_Sendrecv(&HPHY_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_LAST_H_P
            ,&HPHY_LOCAL(t + k,0, 1),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&UPHY_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_LAST_U_P
            ,&UPHY_LOCAL(t + k,0, 1),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_FIRST_U_P, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&VPHY_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_LAST_V_P
            ,&VPHY_LOCAL(t + k,0, 1),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_FIRST_V_P, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&UFIL_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_LAST_U_F
            ,&UFIL_LOCAL(t + k,0, 1),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_FIRST_U_F, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&VFIL_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_LAST_V_F
            ,&VFIL_LOCAL(t + k,0, 1),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_FIRST_V_F, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&HFIL_LOCAL(t + k,1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank-NbCol,TAG_LAST_H_F
            ,&HFIL_LOCAL(t + k,0, 1),size_y/NbCol, MPI_DOUBLE,my_rank-NbCol,TAG_FIRST_H_F, MPI_COMM_WORLD,&status);

            //printf("P#%d:mpirettype_1%d\n",my_rank, mpi_ret_type);

        }
        if(my_rank<NbCol*(NbLi-1)) //envoie celui du bas
        {
            mpi_ret_type = MPI_Sendrecv(&HPHY_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank+NbCol,TAG_FIRST_H_P
            ,&HPHY_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank+NbCol,TAG_LAST_H_P, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&UPHY_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank+NbCol,TAG_FIRST_U_P
            ,&UPHY_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank+NbCol,TAG_LAST_U_P, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&VPHY_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank+NbCol,TAG_FIRST_V_P
            ,&VPHY_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank+NbCol,TAG_LAST_V_P, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&UFIL_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank+NbCol,TAG_FIRST_U_F
            ,&UFIL_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank+NbCol,TAG_LAST_U_F, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&VFIL_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank+NbCol,TAG_FIRST_V_F
            ,&VFIL_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank+NbCol,TAG_LAST_V_F, MPI_COMM_WORLD,&status);

            mpi_ret_type = MPI_Sendrecv(&HFIL_LOCAL(t + k,local_size_x-2, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE, my_rank+NbCol,TAG_FIRST_H_F
            ,&HFIL_LOCAL(t + k,local_size_x-1, my_rank%NbCol!=0),size_y/NbCol, MPI_DOUBLE,my_rank+NbCol,TAG_LAST_H_F, MPI_COMM_WORLD,&status);

            //printf("P#%d:mpirettype_2%d\n",my_rank, mpi_ret_type);
        }

        /*a droite a gauche */
        if(my_rank%NbCol!=0) //tout les processus sauf ceux qui sont sur la colonne de gauche, on envoie la colonne tout a gauche
        {
         	for(i=0;i<size_x/NbLi;i++){
                    MPI_Sendrecv(&HPHY_LOCAL(t + k,i+(my_rank>=NbCol), 1),1, MPI_DOUBLE, my_rank-1,TAG_LAST_H_P
                    ,&HPHY_LOCAL(t + k,i+(my_rank>=NbCol), 0),1, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                    MPI_Sendrecv(&UPHY_LOCAL(t + k,i+(my_rank>=NbCol), 1),1, MPI_DOUBLE, my_rank-1,TAG_LAST_H_P
                    ,&UPHY_LOCAL(t + k,i+(my_rank>=NbCol), 0),1, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
           

                    MPI_Sendrecv(&VPHY_LOCAL(t + k,i+(my_rank>=NbCol), 1),1, MPI_DOUBLE, my_rank-1,TAG_LAST_H_P
                    ,&VPHY_LOCAL(t + k,i+(my_rank>=NbCol), 0),1, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);


                    MPI_Sendrecv(&HFIL_LOCAL(t + k,i+(my_rank>=NbCol), 1),1, MPI_DOUBLE, my_rank-1,TAG_LAST_H_P
                    ,&HFIL_LOCAL(t + k,i+(my_rank>=NbCol), 0),1, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);


                    MPI_Sendrecv(&UFIL_LOCAL(t + k,i+(my_rank>=NbCol), 1),1, MPI_DOUBLE, my_rank-1,TAG_LAST_H_P
                    ,&UFIL_LOCAL(t + k,i+(my_rank>=NbCol), 0),1, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);


                    MPI_Sendrecv(&VFIL_LOCAL(t + k,i+(my_rank>=NbCol), 1),1, MPI_DOUBLE, my_rank-1,TAG_LAST_H_P
                    ,&VFIL_LOCAL(t + k,i+(my_rank>=NbCol), 0),1, MPI_DOUBLE,my_rank-1,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);


            }

            /*petit carre en haut a gauche*/
            if(my_rank>=NbCol) //on envoie celui du haut
            {
                MPI_Sendrecv(&HPHY_LOCAL(t + k,1, 1),1, MPI_DOUBLE, my_rank-1-NbCol,TAG_LAST_H_P
                ,&HPHY_LOCAL(t + k,0, 0),1, MPI_DOUBLE,my_rank-1-NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                MPI_Sendrecv(&UPHY_LOCAL(t + k,1, 1),1, MPI_DOUBLE, my_rank-1-NbCol,TAG_LAST_H_P
                ,&UPHY_LOCAL(t + k,0, 0),1, MPI_DOUBLE,my_rank-1-NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                MPI_Sendrecv(&VPHY_LOCAL(t + k,1, 1),1, MPI_DOUBLE, my_rank-1-NbCol,TAG_LAST_H_P
                ,&VPHY_LOCAL(t + k,0,0),1, MPI_DOUBLE,my_rank-1-NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                MPI_Sendrecv(&HFIL_LOCAL(t + k,1, 1),1, MPI_DOUBLE, my_rank-1-NbCol,TAG_LAST_H_P
                ,&HFIL_LOCAL(t + k,0,0),1, MPI_DOUBLE,my_rank-1-NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                MPI_Sendrecv(&UFIL_LOCAL(t + k,1, 1),1, MPI_DOUBLE, my_rank-1-NbCol,TAG_LAST_H_P
                ,&UFIL_LOCAL(t + k,0,0),1, MPI_DOUBLE,my_rank-1-NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                MPI_Sendrecv(&VFIL_LOCAL(t + k,1, 1),1, MPI_DOUBLE, my_rank-1-NbCol,TAG_LAST_H_P
                ,&VFIL_LOCAL(t + k,0,0),1, MPI_DOUBLE,my_rank-1-NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            }        

            /*petit carre en bas a gauche */
            if(my_rank<NbCol*(NbLi-1))
            {
                MPI_Sendrecv(&HPHY_LOCAL(t + k,1, local_size_y-2),1, MPI_DOUBLE, my_rank-1+NbCol,TAG_LAST_H_P
                ,&HPHY_LOCAL(t + k,0, local_size_y-1),1, MPI_DOUBLE,my_rank-1+NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                MPI_Sendrecv(&UPHY_LOCAL(t + k,1, local_size_y-2),1, MPI_DOUBLE, my_rank-1+NbCol,TAG_LAST_H_P
                ,&UPHY_LOCAL(t + k,0, local_size_y-1),1, MPI_DOUBLE,my_rank-1+NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                MPI_Sendrecv(&VPHY_LOCAL(t + k,1, local_size_y-2),1, MPI_DOUBLE, my_rank-1+NbCol,TAG_LAST_H_P
                ,&VPHY_LOCAL(t + k,0, local_size_y-1),1, MPI_DOUBLE,my_rank-1+NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                MPI_Sendrecv(&HFIL_LOCAL(t + k,1, local_size_y-2),1, MPI_DOUBLE, my_rank-1+NbCol,TAG_LAST_H_P
                ,&HFIL_LOCAL(t + k,0, local_size_y-1),1, MPI_DOUBLE,my_rank-1+NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                MPI_Sendrecv(&UFIL_LOCAL(t + k,1, local_size_y-2),1, MPI_DOUBLE, my_rank-1+NbCol,TAG_LAST_H_P
                ,&UFIL_LOCAL(t + k,0, local_size_y-1),1, MPI_DOUBLE,my_rank-1+NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                MPI_Sendrecv(&VFIL_LOCAL(t + k,1, local_size_y-2),1, MPI_DOUBLE, my_rank-1+NbCol,TAG_LAST_H_P
                ,&VFIL_LOCAL(t + k,0, local_size_y-1),1, MPI_DOUBLE,my_rank-1+NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            }

            //printf("P#%d:mpirettype_1%d\n",my_rank, mpi_ret_type);

        }
        if((my_rank+1)%NbCol!=0) //tout les processus sauf ceux sur la colonne de droite, on envoiel aligne de droite
        {
          	for(i=0;i<size_x/NbLi;i++){
                    MPI_Sendrecv(&HPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2),1, MPI_DOUBLE, my_rank+1,TAG_LAST_H_P
                    ,&HPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1),1, MPI_DOUBLE,my_rank+1,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);

                    MPI_Sendrecv(&UPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2),1, MPI_DOUBLE, my_rank+1,TAG_LAST_H_P
                    ,&UPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1),1, MPI_DOUBLE,my_rank+1,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);

                    MPI_Sendrecv(&VPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2),1, MPI_DOUBLE, my_rank+1,TAG_LAST_H_P
                    ,&VPHY_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1),1, MPI_DOUBLE,my_rank+1,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);

                    MPI_Sendrecv(&HFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2),1, MPI_DOUBLE, my_rank+1,TAG_LAST_H_P
                    ,&HFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1),1, MPI_DOUBLE,my_rank+1,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);

                    MPI_Sendrecv(&UFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2),1, MPI_DOUBLE, my_rank+1,TAG_LAST_H_P
                    ,&UFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1),1, MPI_DOUBLE,my_rank+1,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);

                    MPI_Sendrecv(&VFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-2),1, MPI_DOUBLE, my_rank+1,TAG_LAST_H_P
                    ,&VFIL_LOCAL(t + k,i+(my_rank>=NbCol), local_size_y-1),1, MPI_DOUBLE,my_rank+1,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);

            }

            /*petit carre en haut a droite*/
            if(my_rank>=NbCol) //on envoie celui du haut
            {
                MPI_Sendrecv(&HPHY_LOCAL(t + k,1, local_size_y-2),1, MPI_DOUBLE, my_rank+1-NbCol,TAG_LAST_H_P
                ,&HPHY_LOCAL(t + k,0, local_size_x-1),1, MPI_DOUBLE,my_rank+1-NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
           
                MPI_Sendrecv(&UPHY_LOCAL(t + k,1, local_size_y-2),1, MPI_DOUBLE, my_rank+1-NbCol,TAG_LAST_H_P
                ,&UPHY_LOCAL(t + k,0, local_size_x-1),1, MPI_DOUBLE,my_rank+1-NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(&VPHY_LOCAL(t + k,1, local_size_y-2),1, MPI_DOUBLE, my_rank+1-NbCol,TAG_LAST_H_P
                ,&VPHY_LOCAL(t + k,0, local_size_x-1),1, MPI_DOUBLE,my_rank+1-NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(&HFIL_LOCAL(t + k,1, local_size_y-2),1, MPI_DOUBLE, my_rank+1-NbCol,TAG_LAST_H_P
                ,&HFIL_LOCAL(t + k,0, local_size_x-1),1, MPI_DOUBLE,my_rank+1-NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(&UFIL_LOCAL(t + k,1, local_size_y-2),1, MPI_DOUBLE, my_rank+1-NbCol,TAG_LAST_H_P
                ,&UFIL_LOCAL(t + k,0, local_size_x-1),1, MPI_DOUBLE,my_rank+1-NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);

                MPI_Sendrecv(&VFIL_LOCAL(t + k,1, local_size_y-2),1, MPI_DOUBLE, my_rank+1-NbCol,TAG_LAST_H_P
                ,&VFIL_LOCAL(t + k,0, local_size_x-1),1, MPI_DOUBLE,my_rank+1-NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            }        

            /*petit carre en bas a droite */
            if(my_rank<NbCol*(NbLi-1))
            {
                MPI_Sendrecv(&HPHY_LOCAL(t + k,local_size_x-2, local_size_y-2),1, MPI_DOUBLE, my_rank+1+NbCol,TAG_LAST_H_P
                ,&HPHY_LOCAL(t + k,local_size_x-1, local_size_y-1),1, MPI_DOUBLE,my_rank+1+NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                MPI_Sendrecv(&UPHY_LOCAL(t + k,local_size_x-2, local_size_y-2),1, MPI_DOUBLE, my_rank+1+NbCol,TAG_LAST_H_P
                ,&UPHY_LOCAL(t + k,local_size_x-1, local_size_y-1),1, MPI_DOUBLE,my_rank+1+NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                MPI_Sendrecv(&VPHY_LOCAL(t + k,local_size_x-2, local_size_y-2),1, MPI_DOUBLE, my_rank+1+NbCol,TAG_LAST_H_P
                ,&VPHY_LOCAL(t + k,local_size_x-1, local_size_y-1),1, MPI_DOUBLE,my_rank+1+NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                MPI_Sendrecv(&HFIL_LOCAL(t + k,local_size_x-2, local_size_y-2),1, MPI_DOUBLE, my_rank+1+NbCol,TAG_LAST_H_P
                ,&HFIL_LOCAL(t + k,local_size_x-1, local_size_y-1),1, MPI_DOUBLE,my_rank+1+NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                MPI_Sendrecv(&UFIL_LOCAL(t + k,local_size_x-2, local_size_y-2),1, MPI_DOUBLE, my_rank+1+NbCol,TAG_LAST_H_P
                ,&UFIL_LOCAL(t + k,local_size_x-1, local_size_y-1),1, MPI_DOUBLE,my_rank+1+NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            
                MPI_Sendrecv(&VFIL_LOCAL(t + k,local_size_x-2, local_size_y-2),1, MPI_DOUBLE, my_rank+1+NbCol,TAG_LAST_H_P
                ,&VFIL_LOCAL(t + k,local_size_x-1, local_size_y-1),1, MPI_DOUBLE,my_rank+1+NbCol,TAG_FIRST_H_P, MPI_COMM_WORLD,&status);
            }  
        
        
        }
    }

    mpi_ret_type++;
    //printf("P#%d:line%d\n",my_rank,189);
    for (int j = 0; j < local_size_y; j++) {
      for (int i = 0; i < local_size_x+1; i++) {
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
            HPHY_LOCAL(t, i, j) = hPhy_forward(t, i, j);
            UPHY_LOCAL(t, i, j) = uPhy_forward(t, i, j);
            VPHY_LOCAL(t, i, j) = vPhy_forward(t, i, j);
            HFIL_LOCAL(t, i, j) = hFil_forward(t, i, j);
            UFIL_LOCAL(t, i, j) = uFil_forward(t, i, j);
            VFIL_LOCAL(t, i, j) = vFil_forward(t, i, j);
          }
      }
    }
    //for(k=0;k<2;k++)
    {
        printf("---------------------------- Magic The Gathering ----------------------------\n");
        for(i=0;i<size_x/NbLi;i++)
        {
            MPI_Gather(&HFIL_LOCAL(t,i+(my_rank>=NbCol), (my_rank%NbCol!=0))/*+size_y*(my_rank!=0)*/,size_y/NbCol/*(local_size_x-1-1*(my_rank!=0 && my_rank!=NP-1))*/
            ,MPI_DOUBLE,&HFIL(t, (my_rank%NbCol)*size_y/NbCol, (my_rank%NbLi)*size_x/NbLi),size_y/NbCol/*(local_size_x-1-1*(my_rank!=0 && my_rank!=NP-1))*/,MPI_DOUBLE,0,MPI_COMM_WORLD);


        }


        printf("---------------------------- End of The Gathering ----------------------------\n");
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
  }
  if(my_rank==0)
  if (file_export) {
  	printf("finalize_export\n");
    finalize_export(file);
  }
}

int main(int argc, char **argv)  {

  /* Variables liees au chronometrage */
  double debut=0, fin=0;
  root = 0;
  parse_args(argc, argv);
  printf("Command line options parsed\n");

  alloc();
  printf("Memory allocated\n");
  printf("State initialised\n");
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&NP);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

  /*NP doit etre une puissance de 2*/
  NP_temp = NP;
  local_size_x = size_x;
  local_size_y = size_y;
  NbCol = 0;
  NbLi = 0;
  ligne_colonne = 0; //0 pour ligne 1 pour colonne
  while(NP_temp!=1)
  {
    NP_temp = NP_temp/2;
    if(ligne_colonne==0) //ligne
    {
      local_size_x = local_size_x/2;
      NbLi++;
    }
    else
    {
      local_size_y = local_size_y/2;
      NbCol++;
    }
    ligne_colonne = (ligne_colonne+1)%2;
  }

  if(my_rank>=NbCol)
    local_size_x++;
  if(my_rank<(NbLi-1)*NbCol)
    local_size_x++;
  if(my_rank%NbCol!=0)
    local_size_y++;
  if((my_rank+1)%NbCol!=0)
    local_size_y++;


  gauss_init();
  printf("P#%d:local size x:%d , y:%d\n",my_rank,local_size_x,size_y);
  printf("P#%d:size x:%d , y:%d\n",my_rank,size_x,size_y);

  if(my_rank==0)
  {
    /* debut du chronometrage */
    debut = my_gettimeofday();
    printf("NP:%d\n",NP);
  }

  if(size_x%NP!=0)
  {
    fprintf(stderr,"ERROR:NP ne divise pas size_y\n");
    return EXIT_FAILURE;
  }
  forward_bloc(); 
  printf("State computed\n");

  if(my_rank==0)
  {
     /* fin du chronometrage */
    fin = my_gettimeofday();
    printf("#%d-Temps total de calcul : %g seconde(s) \n",my_rank,fin - debut);
  }

  dealloc();
  printf("Memory freed\n");

  MPI_Finalize();

  return EXIT_SUCCESS;



}