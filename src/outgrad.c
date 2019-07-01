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
 * output of gradients to GRAD_FILE:
 * S. Dunkl 2013 
 -------------------------------------------------------------*/


#include "fd.h"

void outgrad(int nx,int ny,int nz,float ***grad1, float ***grad2,float ***grad3, float finv, int iteration, char outfile[STRING_SIZE]){

	extern int POS[4];
	
	/*extern char  GRAD_FILE[STRING_SIZE];*/
	FILE *fpmod1, *fpmod2, *fpmod3;
	extern int MYID;
	/*extern FILE *FP;*/
	
	int i,j,k;
	char gradfile1[STRING_SIZE],gradfile2[STRING_SIZE],gradfile3[STRING_SIZE];
	char gradfile4[STRING_SIZE],gradfile5[STRING_SIZE],gradfile6[STRING_SIZE];
	
	
		
	sprintf(gradfile1,"%s.vp_%4.2fHz_it%d.%i.%i.%i",outfile,finv,iteration,POS[1],POS[2],POS[3]);
	sprintf(gradfile2,"%s.vs_%4.2fHz_it%d.%i.%i.%i",outfile,finv,iteration,POS[1],POS[2],POS[3]);
	sprintf(gradfile3,"%s.rho_%4.2fHz_it%d.%i.%i.%i",outfile,finv,iteration,POS[1],POS[2],POS[3]);
	
	
	fpmod1=fopen(gradfile1,"w");
	fpmod2=fopen(gradfile2,"w");
	fpmod3=fopen(gradfile3,"w");
	
	for (k=1;k<=nz;k++){
		for (i=1;i<=nx;i++){
			for (j=1;j<=ny;j++){

			fwrite(&grad1[j][i][k], sizeof(float), 1,fpmod1);
			fwrite(&grad2[j][i][k], sizeof(float), 1,fpmod2);
			fwrite(&grad3[j][i][k], sizeof(float), 1,fpmod3);
	
			}
		}
	}

	fclose(fpmod1);
	fclose(fpmod2);
	fclose(fpmod3);
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(MYID==0){
		sprintf(gradfile4,"%s.vp_%4.2fHz_it%d",outfile,finv,iteration);
		mergemod(gradfile4,3);
		sprintf(gradfile5,"%s.vs_%4.2fHz_it%d",outfile,finv,iteration);
		mergemod(gradfile5,3);
		sprintf(gradfile6,"%s.rho_%4.2fHz_it%d",outfile,finv,iteration);
		mergemod(gradfile6,3);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	remove(gradfile1);
	remove(gradfile2);
	remove(gradfile3);
}

	

