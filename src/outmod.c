/*------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS3D.
 * 
 * IFOS3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS3D. See file COPYING and/or 
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/

/*------------------------------------------------------------
 * output of model parameters vp,vs and rho to MOD_OUT_FILE;
 * S. Dunkl 2013 
 -------------------------------------------------------------*/



#include "fd.h"

void outmod(int nx,int ny,int nz,float ***rho, float ***pi,float ***u, int iteration){

	extern int POS[4];
	extern char  MOD_OUT_FILE[STRING_SIZE];
	FILE *fpmod1, *fpmod2, *fpmod3;
	extern int MYID;

	int i,j,k;
	char modfile1[STRING_SIZE],modfile2[STRING_SIZE],modfile3[STRING_SIZE], modfile4[STRING_SIZE];
	char modfile5[STRING_SIZE],modfile6[STRING_SIZE];
	float vp,vs;

	sprintf(modfile1,"%s.vp_it%d.%i.%i.%i",MOD_OUT_FILE,iteration,POS[1],POS[2],POS[3]);
	sprintf(modfile2,"%s.vs_it%d.%i.%i.%i",MOD_OUT_FILE,iteration,POS[1],POS[2],POS[3]);
	sprintf(modfile3,"%s.rho_it%d.%i.%i.%i",MOD_OUT_FILE,iteration,POS[1],POS[2],POS[3]);
	
	fpmod1=fopen(modfile1,"w");
	fpmod2=fopen(modfile2,"w");
	fpmod3=fopen(modfile3,"w");
	
	
	for (k=1;k<=nz;k++){
		for (i=1;i<=nx;i++){	
			for (j=1;j<=ny;j++){

			vp=sqrt(pi[j][i][k]/rho[j][i][k]);
			vs=sqrt(u[j][i][k]/rho[j][i][k]);
			
			  
			fwrite(&vp, sizeof(float), 1,fpmod1);
			fwrite(&vs, sizeof(float), 1,fpmod2);
			fwrite(&rho[j][i][k], sizeof(float), 1,fpmod3);
	
			}
		}
	}
	
	

	fclose(fpmod1);
	fclose(fpmod2);
	fclose(fpmod3);
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(MYID==0){
		sprintf(modfile4,"%s.vp_it%d",MOD_OUT_FILE,iteration);
		mergemod(modfile4,3);
		sprintf(modfile5,"%s.vs_it%d",MOD_OUT_FILE,iteration);
		mergemod(modfile5,3);
		sprintf(modfile6,"%s.rho_it%d",MOD_OUT_FILE,iteration);
		mergemod(modfile6,3);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	remove(modfile1);
	remove(modfile2);
	remove(modfile3);	
}