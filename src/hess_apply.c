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
 * Preconditiioning by applying approximate Hessian to gradient.
 * S. Butzer 2013
 -------------------------------------------------------------------------*/


#include "fd.h"

void hess_apply(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, float ***grad1, float ***grad2,float ***grad3, float ***hess1, float ***hess2, float ***hess3, float finv, int iteration){

	extern float WATER_HESS[3], VP0, VS0, RHO0;
	extern int NXG, NYG, NZG, FW, NX, NY, NZ;
	extern FILE *FP;
	extern int LBFGS;
	
	float wl[3],buf[3];
	int i,j,k,buf1;
	
		
	buf[0]=0.0;buf[1]=0.0;buf[2]=0.0;
	wl[0]=0.0;wl[1]=0.0;wl[2]=0.0;
	
	/*Normalize gradient for LBFGS*/
	if(LBFGS){
		for (j=1;j<=ny2;j++){
			for (i=1;i<=nx2;i++){
				for (k=1;k<=nz2;k++){
				grad1[j][i][k]=grad1[j][i][k]*VP0;
				grad2[j][i][k]=grad2[j][i][k]*VS0;
				grad3[j][i][k]=grad3[j][i][k]*RHO0;
				/*hess1[j][i][k]=hess1[j][i][k]*pow(vp0,2.0);
				hess2[j][i][k]=hess2[j][i][k]*pow(vs0,2.0);
				hess3[j][i][k]=hess3[j][i][k]*pow(rho0,2.0);*/
				}
			}
		}
	}
	
	
	for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){
			for (k=nz1;k<=nz2;k++){
			  wl[0]+=log10(hess1[j][i][k]);
			  wl[1]+=log10(hess2[j][i][k]);
			  wl[2]+=log10(hess3[j][i][k]);
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&wl,&buf,3,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	
	buf1=(NXG-2*FW)*(NYG-2*FW)*(NZG-2*FW);
	wl[0]=1.0*pow(10,(buf[0]/buf1));
	wl[1]=1.0*pow(10,(buf[1]/buf1));
	wl[2]=1.0*pow(10,(buf[2]/buf1));	
	fprintf(FP,"w1=%e, w2=%e, w3=%e \n", wl[0], wl[1], wl[2]);
	
	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
			for (k=1;k<=NZ;k++){
			  grad1[j][i][k]=grad1[j][i][k]/(hess1[j][i][k]+WATER_HESS[0]);
			  grad2[j][i][k]=grad2[j][i][k]/(hess2[j][i][k]+WATER_HESS[1]);
			  grad3[j][i][k]=grad3[j][i][k]/(hess3[j][i][k]+WATER_HESS[2]);
			}
		}
	}
}
