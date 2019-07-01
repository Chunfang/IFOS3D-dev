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

/*------------------------------------------------------------------------
 *   Read elastic model properties (vp,vs,density) from files  
 *  ----------------------------------------------------------------------*/


#include "fd.h"

void readmod(float  ***  rho, float ***  pi, float ***  u, 
float ***  taus, float ***  taup, float *  eta){

	
	extern float DT, *FL, TAU;
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L;
	extern char  MFILE[STRING_SIZE];
	extern FILE *FP;

		
	/* local variables */
	float rhov, muv, piv, vp, vs; /*qp, qs;*/
	float *pts, ts, tp, sumu, sumpi, ws;
	int i, j, l, k, ii, jj, kk;
	FILE *fp_vs, *fp_vp, *fp_rho;/**fp_qp ,*fp_qs;*/
	char filename[STRING_SIZE];

	/* choose data format: ASCII: format=2, BINARY: format=3*/
	const int format=3; 
	

	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}

	ts=TAU;  
	tp=TAU;

	ws=2.0*PI*FL[1];
	
	sumu=0.0; 
	sumpi=0.0;
	for (l=1;l<=L;l++){
		sumu=sumu+((ws*ws*pts[l]*pts[l]*ts)/(1.0+ws*ws*pts[l]*pts[l]));
		sumpi=sumpi+((ws*ws*pts[l]*pts[l]*tp)/(1.0+ws*ws*pts[l]*pts[l]));
	}

		

	   fprintf(FP,"\n...reading model information from model-files...\n");

	   fprintf(FP,"\t P-wave velocities:\n\t %s.vp\n\n",MFILE);
	   sprintf(filename,"%s.vp",MFILE);
	   fp_vp=fopen(filename,"r");
	   if (fp_vp==NULL) err(" Could not open model file for P velocities ! ");


	   fprintf(FP,"\t Shear wave velocities:\n\t %s.vs\n\n",MFILE);
	   sprintf(filename,"%s.vs",MFILE);
	   fp_vs=fopen(filename,"r");
	   if (fp_vs==NULL) err(" Could not open model file for shear velocities ! ");

	   fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
	   sprintf(filename,"%s.rho",MFILE);
	   fp_rho=fopen(filename,"r");
	   if (fp_rho==NULL) err(" Could not open model file for densities ! ");

/*	   fprintf(FP,"\t Qp:\n\t %s.qp\n\n",MFILE);
	   sprintf(filename,"%s.qp",MFILE);
	   fp_qp=fopen(filename,"r");
	   if (fp_qp==NULL) err(" Could not open model file for Qp-values ! ");

	   fprintf(FP,"\t Qs:\n\t %s.qs\n\n",MFILE);
	   sprintf(filename,"%s.qs",MFILE);
	   fp_qs=fopen(filename,"r");
	   if (fp_qs==NULL) err(" Could not open model file for Qs-values ! ");
*/	   

	/* loop over global grid */
		for (k=1;k<=NZG;k++){ 
			for (i=1;i<=NXG;i++){
				for (j=1;j<=NYG;j++){

				vp=readdsk(fp_vp, format);
				vs=readdsk(fp_vs, format);
				rhov=readdsk(fp_rho , format);
				/*qp=readdsk(fp_qp, format);
				qs=readdsk(fp_qs, format);*/
				
				
				muv=vs*vs*rhov/(1.0+sumu);
				piv=vp*vp*rhov/(1.0+sumpi);

				/* only the PE which belongs to the current global gridpoint 
				is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY)) && 
				    (POS[3]==((k-1)/NZ))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
					kk=k-POS[3]*NZ;
					
					if(L){
					taus[jj][ii][kk]=ts;
					taup[jj][ii][kk]=tp;}
					
					u[jj][ii][kk]=muv;
					rho[jj][ii][kk]=rhov;
					pi[jj][ii][kk]=piv;
					
			
				}
				}
			}
		}




	fclose(fp_vp);
	fclose(fp_vs);
	fclose(fp_rho);
/*	fclose(fp_qp);
	fclose(fp_qs);*/
	
	
	/* each PE writes his model to disk */
	   
	   
	/*sprintf(filename,"%s.IFOS.pwavemod",MFILE);
	writemod(filename,pwavemod,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(filename,3);
	
	sprintf(filename,"%s.IFOS.swavemod",MFILE);
	writemod(filename,swavemod,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(filename,3); 
	
	sprintf(filename,"%s.IFOS.rho",MFILE);
	writemod(filename,rho,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(filename,3);*/

	free_vector(pts,1,L);
}




