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

/* -----------------------------------------------------------------------
 *  generation of an elastic or viscoelastic model specified by vp, vs and rho
 -------------------------------------------------------------------------*/


#include "fd.h"

void model(float  ***  rho, float ***  pi, float ***  u, 
float ***  taus, float ***  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, *FL, TAU;
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
	extern FILE *FP;
	/* local variables */
	float muv, piv, ws;
	float Vp, Vs, Rho;
	float *pts=NULL, ts=0.0, tp=0.0, sumu, sumpi;
	int i, j, k, l, ii, jj, kk;

	/* parameters for layer 1 */
	const float vp1=6400.0, vs1=4100.0, rho1=2800.0;
	const float vp2=6400.0, vs2=4100.0, rho2=2800.0;
      /*const float vp2=7000.0, vs2=3900.0, rho2=2800.0;*/

	/*-----------------------------------------------------------------------*/
	sumu=0.0; 
	sumpi=0.0;

	fprintf(FP," start model creation");
	
	 if(L){
		/* vector for maxwellbodies */
		pts=vector(1,L);
		for (l=1;l<=L;l++) {
			pts[l]=1.0/(2.0*PI*FL[l]);
			eta[l]=DT/pts[l];
		}
		ts=TAU;  
		tp=TAU;
		ws=2.0*PI*FL[1];
		for (l=1;l<=L;l++){
			sumu=sumu+((ws*ws*pts[l]*pts[l]*ts)/(1.0+ws*ws*pts[l]*pts[l]));
			sumpi=sumpi+((ws*ws*pts[l]*pts[l]*tp)/(1.0+ws*ws*pts[l]*pts[l]));
		}
	 }

	/* loop over global grid */
	for (j=1;j<=NYG;j++){
		for (i=1;i<=NXG;i++){
			for (k=1;k<=NZG;k++){

				
				Vp=vp1; Vs=vs1; Rho=rho1; 

				muv=Vs*Vs*Rho/(1.0+sumu);
				piv=Vp*Vp*Rho/(1.0+sumpi);

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
						taup[jj][ii][kk]=tp;
					}
					u[jj][ii][kk]=muv;
					rho[jj][ii][kk]=Rho;
					pi[jj][ii][kk]=piv;
				}
			}
		}
	}

	for (j=50;j<=80;j++){
		for (i=55;i<=85;i++){
			for (k=50;k<=80;k++){

			/*for (j=(NYG/2+1);j<=NYG;j++){*/

				Vp=0.0; Vs=0.0; Rho=0.0;
				Vp=vp2; Vs=vs2; Rho=rho2; 

				muv=Vs*Vs*Rho/(1.0+sumu);
				piv=Vp*Vp*Rho/(1.0+sumpi);

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
						taup[jj][ii][kk]=tp;
					}
					u[jj][ii][kk]=muv;
					rho[jj][ii][kk]=Rho;
					pi[jj][ii][kk]=piv;
				}
			}
			
		}
	}	


	MPI_Barrier(MPI_COMM_WORLD);

         if(L)	free_vector(pts,1,L);
}



