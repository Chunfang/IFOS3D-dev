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

/*--------------------------------------------------------------------------
 * calculates the conjugate gradient direction (MORA, 1987)
 * S. Butzer 2015
 * ------------------------------------------------------------------------*/

#include "fd.h"

void conjugate(int nx,int ny,int nz, float ***grad1, float ***grad2,float ***grad3, float ***gradprior1, float ***gradprior2, float ***gradprior3, float ***gradprior4, float ***gradprior5, float ***gradprior6,float *beta, int iteration,int cdf){

	extern FILE *FP;

	int i,j,k;
	float betadum[6],betadum1[6];
	
	betadum[0]=0.0; betadum[1]=0.0;
	betadum[2]=0.0; betadum[3]=0.0;
	betadum[4]=0.0; betadum[5]=0.0;
	
	beta[0]=0.0; beta[1]=0.0; beta[2]=0.0;

	fprintf(FP,"\n Calculate conjugate gradient ");
	
	if(cdf==1){
		for (j=1;j<=ny;j++){
			for (i=1;i<=nx;i++){
				for (k=1;k<=nz;k++){
				gradprior1[j][i][k]=grad1[j][i][k]; 
				gradprior2[j][i][k]=grad2[j][i][k];
				gradprior3[j][i][k]=grad3[j][i][k];
				gradprior4[j][i][k]=grad1[j][i][k]; 
				gradprior5[j][i][k]=grad2[j][i][k];
				gradprior6[j][i][k]=grad3[j][i][k];
				}
			}
		}
	}

	else{
		/*Calculation of beta for conjugate gradient*/
		for (j=1;j<=ny;j++){
			for (i=1;i<=nx;i++){
				for (k=1;k<=nz;k++){		
					betadum[0]+=grad1[j][i][k]*(grad1[j][i][k]-gradprior1[j][i][k]);
					betadum[1]+=gradprior1[j][i][k]*gradprior1[j][i][k];
					betadum[2]+=grad2[j][i][k]*(grad2[j][i][k]-gradprior2[j][i][k]);
					betadum[3]+=gradprior2[j][i][k]*gradprior2[j][i][k];
					betadum[4]+=grad3[j][i][k]*(grad3[j][i][k]-gradprior3[j][i][k]);
					betadum[5]+=gradprior3[j][i][k]*gradprior3[j][i][k];
					
					/*save preconditioned gradient for next iteration*/
					gradprior1[j][i][k]=grad1[j][i][k]; 
					gradprior2[j][i][k]=grad2[j][i][k];
					gradprior3[j][i][k]=grad3[j][i][k];			
				}
			}
		}

		MPI_Allreduce(&betadum,&betadum1,6,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		
		if(cdf==0){
		if(fabs(betadum1[1])>0.0) beta[0]=betadum1[0]/betadum1[1];
		if(fabs(betadum1[3])>0.0) beta[1]=betadum1[2]/betadum1[3];
		if(fabs(betadum1[5])>0.0) beta[2]=betadum1[4]/betadum1[5];
		fprintf(FP,"\n beta[1]=%e, beta[2]=%e, beta[3]=%e\n",beta[0],beta[1],beta[2] );}
		
		if(beta[0]<0)beta[0]=0.0;
		if(beta[0]<0)beta[1]=0.0;
		if(beta[0]<0)beta[2]=0.0;
		if(beta[0]>0.5)beta[0]=0.5;
		if(beta[0]>0.5)beta[1]=0.5;
		if(beta[0]>0.5)beta[2]=0.5;
		
		/*Calculation of conjugate gradient direction (not normalised)*/
		for (j=1;j<=ny;j++){
			for (i=1;i<=nx;i++){
				for (k=1;k<=nz;k++){
					grad1[j][i][k]=gradprior4[j][i][k]*beta[0]+grad1[j][i][k];
					grad2[j][i][k]=gradprior5[j][i][k]*beta[1]+grad2[j][i][k];
					grad3[j][i][k]=gradprior6[j][i][k]*beta[2]+grad3[j][i][k];
					
					/*save conjugate gradient for next iteration*/
					gradprior4[j][i][k]=grad1[j][i][k]; 
					gradprior5[j][i][k]=grad2[j][i][k];
					gradprior6[j][i][k]=grad3[j][i][k];
				}
			}
		}
	}
		
}
