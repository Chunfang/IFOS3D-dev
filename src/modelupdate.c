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
 * update model for next iteration
 * S. Butzer 2013
 * ------------------------------------------------------------------------*/

#include "fd.h"
void modelupdate(int nx, int ny, int nz, float ***gradvp, float ***gradvs, float ***gradrho, float  ***  rho, float ***  pi, float ***  u, float **bfgsmod1, float step, float *beta, int it_group){

	extern FILE *FP;
	extern int NX, NY, NZ, LBFGS, BFGSNUM, NUMPAR;
	extern float VP0, VS0, RHO0, WEIGHT[3];
	

	int j,i,k,l,l1;
	float vpnew,vsnew,rhonew;
	float max[3],buf[3],max1[3],dummy[3];
	float vp,vs;
	int w=0;
	float scale1=0.0, scale2=0.0, scale3=0.0;
	
	buf[0]=0.0;buf[1]=0.0;buf[2]=0.0;
	max[0]=0.0;max[1]=0.0;max[2]=0.0;
	dummy[0]=0.0; dummy[1]=0.0; dummy[2]=0.0;
	max1[0]=0.0;max1[1]=0.0;max1[2]=0.0;
	
	scale1=WEIGHT[0];
	scale2=WEIGHT[1];
	scale3=WEIGHT[2];
	/*if(iteration<20)scale1=0.7;
	if(iteration>40)scale2=0.7;
	if(iteration>60)scale2=0.4;*/
	
	if(LBFGS) {w=it_group%BFGSNUM;
		  if(w==0) w=BFGSNUM;}

	 fprintf(FP,"\n Message from modelupdate.c ");
	
	/*find and exchange model maxima*/
	for (j=1;j<=ny;j++){
		for (i=1;i<=nx;i++){
			for (k=1;k<=nz;k++){
			vp=0.0;
			vs=0.0;
			vp=sqrt(pi[j][i][k]/rho[j][i][k]);
			vs=sqrt(u[j][i][k]/rho[j][i][k]);
			
			if (vp>max[0]) max[0]=vp;
			if (vs>max[1]) max[1]=vs;
			if (rho[j][i][k]>max[2]) max[2]=rho[j][i][k];
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&max,&buf,3,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	max[0]=buf[0];
	max[1]=buf[1];
	max[2]=buf[2];
	fprintf(FP,"\n Maximum values in current model:\n");
	fprintf(FP," rhomax=%4.2f, vpmax=%4.2f, vsmax=%4.2f \n", max[2], max[0], max[1]);
	
	/*if(iteration==1){vs0=vs0*vs0;
	vp0=vp0*vp0;}
	 else {vs0=1.0;
	 vp0=1.0;}*/

        if(!LBFGS){
		  /*Normalise gradients to maximum*/
		  for (j=1;j<=ny;j++){
			  for (i=1;i<=nx;i++){
				  for (k=1;k<=nz;k++){
				    
				  if (fabs(gradvp[j][i][k])>max1[0]) max1[0]=fabs(gradvp[j][i][k]);
				  if (fabs(gradvs[j][i][k])>max1[1]) max1[1]=fabs(gradvs[j][i][k]);
				  if (fabs(gradrho[j][i][k])>max1[2]) max1[2]=fabs(gradrho[j][i][k]);
				    
				  }
			  }
		  }
		  MPI_Barrier(MPI_COMM_WORLD);
		  MPI_Allreduce(&max1,&dummy,3,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
		  MPI_Barrier(MPI_COMM_WORLD);
		  
		  max1[0]=dummy[0];
		  max1[1]=dummy[1];
		  max1[2]=dummy[2];
		  if(max1[2]==0.0) max1[2]=1.0;
		  //fprintf(FP," max1[0]=%e,max1[1]=%e  \n",max1[0],max1[1]);
		  
		  if (max1[0]>0.0 && max1[1]>0.0){
			  for (j=1;j<=ny;j++){
				  for (i=1;i<=nx;i++){
					  for (k=1;k<=nz;k++){
						  gradvp[j][i][k]= gradvp[j][i][k]/max1[0];
						  gradvs[j][i][k]= gradvs[j][i][k]/max1[1];
						  gradrho[j][i][k]=gradrho[j][i][k]/max1[2];
					  }
				  }
			  }
		  }
		  else fprintf(FP,"Warning: One of the gradients is zero: max1*max2*max3 =0 \n");
	}
	/*if(it_group==1){scale1=4.0; scale2=4.0;}*/
	l=0;
	l1=NX*NY*NZ;
	for (j=1;j<=ny;j++){
		for (i=1;i<=nx;i++){
			for (k=1;k<=nz;k++){

			
			/*update model*/
			vpnew=0.0;
			/*vpnew=sqrt(pi[j][i][k]/rho[j][i][k])+max[0]*step*gradconvp[j][i][k];*/
			vpnew=sqrt(pi[j][i][k]/rho[j][i][k])+step*scale1*gradvp[j][i][k]*VP0;
			
			vsnew=0.0;
			/*vsnew=sqrt(u[j][i][k]/rho[j][i][k])+max[1]*step*gradconvs[j][i][k];*/
			vsnew=sqrt(u[j][i][k]/rho[j][i][k])+step*scale2*gradvs[j][i][k]*VS0;
		      
			rhonew=0.0;
			rhonew=rho[j][i][k]+scale3*step*gradrho[j][i][k]*RHO0;
			
			    if(LBFGS){
			/*save normalised model differences for LBFGS*/
			 l++; l1++;
			bfgsmod1[w][l]=(vpnew-sqrt(pi[j][i][k]/rho[j][i][k]))/VP0;
			if(NUMPAR>1){
			bfgsmod1[w][l1]=(vsnew-sqrt(u[j][i][k]/rho[j][i][k]))/VS0;
			}
			}
			  
			u[j][i][k]=rhonew*vsnew*vsnew;
			pi[j][i][k]=rhonew*vpnew*vpnew;
			rho[j][i][k]=rhonew;
			
			}
		}
	}
		
}
