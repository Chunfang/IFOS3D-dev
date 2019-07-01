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

/* --------------------------------------------------------------------
   Reading (distributed) source positions, timeshift, centre frequency 
   and amplitude from SOURCE_FILE.
----------------------------------------------------------------------*/

#include "fd.h"

void sources(FILE * fpsrc, int *nsrc, float **srcpos){

	/* declaration of extern variables */
	extern int MYID,SRC_MF,SOURCE_TYPE;
	extern FILE *FP;
	extern float TS,REFSRC[3],SRCTSHIFT,FC,AMP;

	int l;
	float xsrc, ysrc, zsrc, tshift, fc=0.0;
	char cline[256];

	if (MYID==0){
			/*srcpos=fmatrix(1,6,1,*nsrc);*/
			/* stype=(int *)malloc(*nsrc*sizeof(int)); */
			
			for (l=1;l<=*nsrc;l++){
				fgets(cline,255,fpsrc);
				switch(sscanf(cline,"%f%f%f%f%f%f%f",&xsrc, &ysrc, &zsrc, &tshift, &srcpos[5][l], &srcpos[6][l], &srcpos[7][l])){
					case 0: xsrc=0.0;
					case 1: zsrc=0.0;
					case 2: ysrc=0.0;
					case 3: tshift=0.0;
					case 4: srcpos[5][l]=FC;
					case 5: srcpos[6][l]=AMP;
					case 6: srcpos[7][l]=SOURCE_TYPE;
				}
				/*change stype to srcpos[7]*/
				/*"y" is used for the vertical coordinate*/
				if(SRC_MF==1) { /* feet */
					srcpos[1][l]=xsrc/0.3048-REFSRC[0];
					srcpos[2][l]=ysrc/0.3048-REFSRC[2];
					srcpos[3][l]=zsrc/0.3048-REFSRC[1];
				}
				else {
					srcpos[1][l]=xsrc-REFSRC[0];
					srcpos[2][l]=ysrc-REFSRC[2];
					srcpos[3][l]=zsrc-REFSRC[1];
				}
				srcpos[4][l]=tshift+SRCTSHIFT*(l-1);
			}
			fclose(fpsrc);

			/* Compute maximum frequency */
			for (l=1;l<=*nsrc;l++)
				if (srcpos[5][l]>fc) fc=srcpos[5][l];
			fprintf(FP," Maximum frequency defined in source file: %6.2f Hz\n",fc);
			TS=1.0/fc;

	}

	/*if (MYID!=0) srcpos=fmatrix(1,6,1,*nsrc);*/
	
	MPI_Bcast(&srcpos[1][1],(*nsrc)*7,MPI_FLOAT,0,MPI_COMM_WORLD);
	
	/*if (MYID!=0) stype=(int *)malloc(*nsrc*sizeof(int)); */
	
	/*MPI_Bcast(&stype[1],*nsrc,MPI_INT,0,MPI_COMM_WORLD);*/
	
	if (MYID==0){
		printf("\n **Message from function source (written by PE %d):\n",MYID);
		printf(" Number of global source positions found: %i\n",*nsrc);
		printf(" x\t\ty\t\tz\t\ttshift\t\tfc\t\tamp\t\tstype\n");}

        MPI_Barrier(MPI_COMM_WORLD);

	/* outputs all sources per each subdomain / node*/
		
	if (MYID==0){
	   for (l=1;l<=*nsrc;l++)
	      printf(" %6.2e\t%6.2e\t%6.2e\t%6.2e\t%6.2e\t%6.2e\t%6.2f\n",
		srcpos[1][l],srcpos[2][l],srcpos[3][l],srcpos[4][l],srcpos[5][l],srcpos[6][l],srcpos[7][l]);
	      printf("\n\n");
	}
        MPI_Barrier(MPI_COMM_WORLD);

}

void pwsources(int *nsrc, float **srcpos){ /* plane wave excitation */

	/* declaration of extern variables */
	extern float PLANE_WAVE_DEPTH, TS, DX, DY, PHI;
	extern int MYID, NXG, NYG, SRCREC, FW,SOURCE_TYPE;
	extern FILE *FP;	

	float x, y, z, tan_phi;
	int  k, l, isrc=0, ixend, iyend;

	if (MYID==0){
		
		if (SRCREC){ /* if SRCREC=1 -> read source positions from file */
			fprintf(FP,"\n Source file is ignored: Plane wave excitation, only. \n");
		} 
		fprintf(FP," Computing source nodes for plane wave excitation.\n");
		fprintf(FP," depth= %5.2f meter, incidence angle= %5.2f degrees.\n",PLANE_WAVE_DEPTH, PHI);

		tan_phi=tan(PHI*PI/180.0);
		fprintf(FP," Message from function sources (written by PE %d):\n",MYID);				
		
	
		/*code relic, not quite sure about the purpose */
		/*if (PHI==0.0) ifw=1;
		if (PHI==0.0) ixend=NXG;
		else ixend=iround((((float)(FW*DX))+(((float)NYG*DY-((float)(FW*DX))-PLANE_WAVE_DEPTH)/tan_phi))/DX);
		
		
		if (ixend>(NXG-ifw+1)) ixend=NXG-ifw+1; */
		
		ixend=NXG-FW;
		iyend=NYG-FW;
		
		/*srcpos=fmatrix(1,6,1,*nsrc);*/
		
		
		/*read from ifos.c */
		/*fprintf(FP," Number of source positions: %i\n",*nsrc);*/
		fprintf(FP," x-range for plane wave: %d to %d gridpoints. \n",FW,ixend);
		fprintf(FP," y-range for plane wave: %d to %d gridpoints. \n",FW,iyend);
						
		/*srcpos=fmatrix(1,6,1,*nsrc);*/
		isrc=0;
		
		for (k=FW;k<=iyend;k++) for (l=FW;l<=ixend;l++){
			x=(float)l*DX;
			y=(float)k*DY;
			z=PLANE_WAVE_DEPTH+(tan_phi*x);
			isrc++;
			srcpos[1][isrc]=x;
			srcpos[2][isrc]=y;
			srcpos[3][isrc]=z;
			srcpos[4][isrc]=0.0;
			srcpos[5][isrc]=1.0/TS;
			srcpos[6][isrc]=1.0;
			srcpos[7][isrc]=SOURCE_TYPE;
		}
		
		/*double check if number of receivers match, counted within loop*/
		fprintf(FP," Number of source positions: %i\n",isrc);
	}
	
	
	/*if (MYID!=0) srcpos=fmatrix(1,6,1,*nsrc);*/
	
	/*MPI_Barrier(MPI_COMM_WORLD);	*/
	
	MPI_Bcast(&srcpos[1][1],(*nsrc)*7,MPI_FLOAT,0,MPI_COMM_WORLD);	
	/*MPI_Bcast(&stype[1],*nsrc,MPI_INT,0,MPI_COMM_WORLD);*/
        MPI_Barrier(MPI_COMM_WORLD);	

		
	/*if (MYID==0){
		fprintf(FP,"\n **Message from function source (written by PE %d):\n",MYID);
		fprintf(FP," Number of global source positions found: %i\n",*nsrc);
		fprintf(FP," MYID\tx\t\ty\t\tz\t\ttshift\t\tfc\t\tamp\n");
		for (l=1;l<=*nsrc;l++)
			fprintf(FP," %6.2e\t%6.2e\t%6.2e\t%6.2e\t%6.2e\t%6.2e\n",
					srcpos[1][l],srcpos[2][l],srcpos[3][l],srcpos[4][l],srcpos[5][l],srcpos[6][l]);
		fprintf(FP,"\n\n");
	}*/

}
