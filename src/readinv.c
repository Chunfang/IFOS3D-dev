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

/*--------------------------------------------------------------------
 *Read inversion parameters from workflow
 --------------------------------------------------------------------*/

#include "fd.h"
void readinv(float *finv, int *nf, int *groupnum,int *itpergroup,int nfmax){

	extern char INV_FILE[STRING_SIZE];
	extern FILE *FP;
	extern int MYID, VERBOSE;
	char buffer[256];
	int mute,i;
	float twin[2],owin[2];
	float abort;
	
	
	FILE *fpinv=NULL;
	
	*groupnum+=1;
	fprintf(FP,"\n\n ********** INVERSION STAGE NR. %d ********** \n",*groupnum);
	if (VERBOSE) fprintf(FP," INV_FILE:%s\n",INV_FILE);
	
	
	
	if(MYID==0){
		fpinv=fopen(INV_FILE,"r");
		if(fpinv==NULL) err("INV_File  could not be opened!");
		
		for(i=1;i<*groupnum;i++){
			fgets(buffer, 255,fpinv);
		}
		
	
	
	fscanf(fpinv,"%d%d%f%d%f%f%f%f%d",&itpergroup[0],&itpergroup[1],&abort,&mute,&twin[0],&twin[1],&owin[0],&owin[1],nf);
	

	fprintf(FP,"\n\n itmin=%d; itmax=%d; \n",itpergroup[0],itpergroup[1]);
	fprintf(FP,"\n Number of Frequencies in Stage Nr.%d: nf=%d \n\n",*groupnum,*nf);
	//fprintf(FP,"\n\n itmin=%d;itmax=%d;abort=%4.2f;mute=%d; twinmin=%4.2f; twinmax=%4.2f;owinmin=%4.2f; owinmax=%4.2f; nf=%d \n\n",itpergroup[0],itpergroup[1],abort,mute,twin[0],twin[1],owin[0],owin[1],*nf);
	
	for(i=0;i<=*nf-1;i++){
		fscanf(fpinv,"%f",&finv[i]);
		fprintf(FP," finv[%d]=%4.2f; ",i,finv[i]);
	}
	
	
	fclose(fpinv);
	i=*nf;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(itpergroup,2,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(nf,1,MPI_INT,0,MPI_COMM_WORLD);
	if(nfmax>1){
	MPI_Bcast(finv,nfmax,MPI_FLOAT,0,MPI_COMM_WORLD);}
	else MPI_Bcast(&finv[0],1,MPI_FLOAT,0,MPI_COMM_WORLD);
	fprintf(FP,"\n");
	/*MPI_Bcast(&abort,1,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&mute,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&owin,2,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&twin,2,MPI_FLOAT,0,MPI_COMM_WORLD);
	*/
	
	
	
		

}