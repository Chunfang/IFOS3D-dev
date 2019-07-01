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

/*---------------------------------------------------------------------------------
 * gradient preconditioning around sources, receivers and at model boundaries:
 * S. Butzer 2013
------------------------------------------------------------------------------------ */


#include "fd.h"

void precon_grad(int nx,int ny,int nz, float ***grad1, float ***grad2,float ***grad3, int nsrc, float **  srcpos, int ntr_glob, int **recpos, float finv, int iteration,int cdf){

	extern float DX, DY, DZ;
	extern int POS[4], NXG, NYG, NZG, FW, DAMPTYPE, FREE_SURF;
	extern FILE *FP;

	int i,j,k,l,n,ii=0,kk=0,jj=0,sx=0,sy=0,sz=0,rx,ry,rz;//,h;
	int sh, rh; /*distance of source and receiver arrays to model edge (grid points)*/
	float r=0.0, G[3];
	G[0]=0.0; G[1]=0.0; G[2]=0.0;
	sh=FW+40;
	rh=FW+40;
	/*h=FW+40;*/
	
	
	fprintf(FP,"\n Gradient preconditioning\n");
	
	if(DAMPTYPE==1){
		for (j=1;j<=ny;j++){
		      for (i=1;i<=nx;i++){
			      for (k=1;k<=nz;k++){
				
				ii=i+POS[1]*nx;
				jj=j+POS[2]*ny;
				kk=k+POS[3]*nz;	 
				
				for (n=1;n<=ntr_glob;n++){
					rx=recpos[1][n];
					ry=recpos[2][n];
					rz=recpos[3][n];
					r=0.0;
					r=sqrt((jj-ry)*(jj-ry)+(ii-rx)*(ii-rx)+(kk-rz)*(kk-rz));
					G[0]=0.0; G[1]=0.0;
					G[0]=1/(1+1000*exp(-0.7*r*r));
					G[1]=1/(1+1000*exp(-0.7*r*r));
					
					grad1[j][i][k]=grad1[j][i][k]*G[0];
					grad2[j][i][k]=grad2[j][i][k]*G[1];
					grad3[j][i][k]=grad3[j][i][k]*G[1];
				}
				/*damping CPML-boundary*/
				/*damping y-direction*/
				G[0]=0.0;
				if(!FREE_SURF){
					r=0.0;
					r=(jj)*1.0;
					G[0]=1/(1+2000*exp(-0.6*r)); 
				}
				
				r=0.0;
				r=(NYG-jj)*1.0;
				G[1]=0.0;
				G[1]=1/(1+2000*exp(-0.6*r));
				
				grad1[j][i][k]=grad1[j][i][k]*G[0]*G[1];
				grad2[j][i][k]=grad2[j][i][k]*G[0]*G[1];
				grad3[j][i][k]=grad3[j][i][k]*G[0]*G[1];
				
				/*damping x-direction*/
				r=0.0;
				r=(ii)*1.0;
				G[0]=0.0;
				G[0]=1/(1+2000*exp(-0.6*r)); 
				
				r=0.0;
				r=(NXG-ii)*1.0;
				G[1]=0.0;
				G[1]=1/(1+2000*exp(-0.6*r));
				
				grad1[j][i][k]=grad1[j][i][k]*G[0]*G[1];
				grad2[j][i][k]=grad2[j][i][k]*G[0]*G[1];
				grad3[j][i][k]=grad3[j][i][k]*G[0]*G[1];
				
				/*damping z-direction*/
				r=0.0;
				r=(kk)*1.0;
				G[0]=0.0;
				G[0]=1/(1+2000*exp(-0.6*r)); 
				
				r=0.0;
				r=(NZG-kk)*1.0;
				G[1]=0.0;
				G[1]=1/(1+2000*exp(-0.6*r));
				
				grad1[j][i][k]=grad1[j][i][k]*G[0]*G[1];
				grad2[j][i][k]=grad2[j][i][k]*G[0]*G[1];
				grad3[j][i][k]=grad3[j][i][k]*G[0]*G[1];
				/*----------------------------------------------*/
				

			      }
		      }
		}  
	}
	
	
	if(DAMPTYPE==2){
	  /* circular damping around source and receiver positions and damping in model boundaries */
		for (j=1;j<=ny;j++){
			for (i=1;i<=nx;i++){
				for (k=1;k<=nz;k++){
				  
				ii=i+POS[1]*nx;
				jj=j+POS[2]*ny;
				kk=k+POS[3]*nz;	  
				
				for (l=1;l<=nsrc;l++) {
					sx=iround(srcpos[1][l]/DX);
					sy=iround(srcpos[2][l]/DY);
					sz=iround(srcpos[3][l]/DZ);
					
					r=0.0;
					r=sqrt((jj-sy)*(jj-sy)+(ii-sx)*(ii-sx)+(kk-sz)*(kk-sz));
					
					G[0]=0.0; G[1]=0.0;		
					G[0]=1/(1+1000*exp(-0.7*r));	  
					G[1]=1/(1+1000*exp(-0.7*r));
					grad1[j][i][k]=grad1[j][i][k]*G[0];
					grad2[j][i][k]=grad2[j][i][k]*G[1];
					grad3[j][i][k]=grad3[j][i][k]*G[1];
				}
				
				for (n=1;n<=ntr_glob;n++){
					rx=recpos[1][n];
					ry=recpos[2][n];
					rz=recpos[3][n];
					r=0.0;
					r=sqrt((jj-ry)*(jj-ry)+(ii-rx)*(ii-rx)+(kk-rz)*(kk-rz));
					G[0]=0.0; G[1]=0.0;
					G[0]=1/(1+1000*exp(-0.7*r*r));
					G[1]=1/(1+1000*exp(-0.7*r*r));
					
					grad1[j][i][k]=grad1[j][i][k]*G[0];
					grad2[j][i][k]=grad2[j][i][k]*G[1];
					grad3[j][i][k]=grad3[j][i][k]*G[1];
				}
				/*damping CPML-boundary*/
				/*damping y-direction*/
				G[0]=0.0;
				if(!FREE_SURF){
					r=0.0;
					r=(jj)*1.0;
					G[0]=1/(1+2000*exp(-0.6*r)); 
				}
				
				r=0.0;
				r=(NYG-jj)*1.0;
				G[1]=0.0;
				G[1]=1/(1+2000*exp(-0.6*r));
				
				grad1[j][i][k]=grad1[j][i][k]*G[0]*G[1];
				grad2[j][i][k]=grad2[j][i][k]*G[0]*G[1];
				grad3[j][i][k]=grad3[j][i][k]*G[0]*G[1];
				
				/*damping x-direction*/
				r=0.0;
				r=(ii)*1.0;
				G[0]=0.0;
				G[0]=1/(1+2000*exp(-0.6*r)); 
				
				r=0.0;
				r=(NXG-ii)*1.0;
				G[1]=0.0;
				G[1]=1/(1+2000*exp(-0.6*r));
				
				grad1[j][i][k]=grad1[j][i][k]*G[0]*G[1];
				grad2[j][i][k]=grad2[j][i][k]*G[0]*G[1];
				grad3[j][i][k]=grad3[j][i][k]*G[0]*G[1];
				
				/*damping z-direction*/
				r=0.0;
				r=(kk)*1.0;
				G[0]=0.0;
				G[0]=1/(1+2000*exp(-0.6*r)); 
				
				r=0.0;
				r=(NZG-kk)*1.0;
				G[1]=0.0;
				G[1]=1/(1+2000*exp(-0.6*r));
				
				grad1[j][i][k]=grad1[j][i][k]*G[0]*G[1];
				grad2[j][i][k]=grad2[j][i][k]*G[0]*G[1];
				grad3[j][i][k]=grad3[j][i][k]*G[0]*G[1];
				/*----------------------------------------------*/
				
				}
			}
		}
	}
	
	if(DAMPTYPE==3){
	  /*damping source and receiver plane and model boundaries (used for transmission geometry)*/
		for (j=1;j<=ny;j++){
			for (i=1;i<=nx;i++){
				for (k=1;k<=nz;k++){		
				/*-------------------------------------------*/
				/*damping y-direction*/
				//h=0;
				//h=iround(srcpos[2][1]/DY);
				r=0.0;
				r=(jj-sh)*1.0;
				if(r<0)	G[0]=0.0;
				else   G[0]=1/(1+1000*exp(-0.7*r*r)); 
				
				r=0.0;
				r=(NYG-jj-rh)*1.0;
				if(r<0)	G[1]=0.0;
				else   G[1]=1/(1+1000*exp(-0.7*r));
				
				grad1[j][i][k]=grad1[j][i][k]*G[0]*G[1];
				grad2[j][i][k]=grad2[j][i][k]*G[0]*G[1];
				grad3[j][i][k]=grad3[j][i][k]*G[0]*G[1];
				
				/*damping x-direction*/
				r=0.0;
				r=(ii)*1.0;
				if(r<0)	G[0]=0.0;
				else   G[0]=1/(1+2000*exp(-0.6*r)); 
				
				r=0.0;
				r=(NXG-ii)*1.0;
				if(r<0)	G[1]=0.0;
				else   G[1]=1/(1+2000*exp(-0.6*r));
				
				grad1[j][i][k]=grad1[j][i][k]*G[0]*G[1];
				grad2[j][i][k]=grad2[j][i][k]*G[0]*G[1];
				grad3[j][i][k]=grad3[j][i][k]*G[0]*G[1];
				
				/*damping z-direction*/
				r=0.0;
				r=(kk)*1.0;
				if(r<0)	G[0]=0.0;
				else   G[0]=1/(1+2000*exp(-0.6*r)); 
				
				r=0.0;
				r=(NZG-kk)*1.0;
				if(r<0)	G[1]=0.0;
				else   G[1]=1/(1+2000*exp(-0.6*r));
				
				grad1[j][i][k]=grad1[j][i][k]*G[0]*G[1];
				grad2[j][i][k]=grad2[j][i][k]*G[0]*G[1];
				grad3[j][i][k]=grad3[j][i][k]*G[0]*G[1];
				/*----------------------------------------------*/
				}
			}
		}
	}
	
	if(DAMPTYPE==0){
	  /*damping only model boundaries*/
		for (j=1;j<=ny;j++){
			for (i=1;i<=nx;i++){
				for (k=1;k<=nz;k++){		
				/*damping CPML-boundary*/
				/*damping y-direction*/
				G[0]=0.0;
				if(!FREE_SURF){
					r=0.0;
					r=(jj)*1.0;
					G[0]=1/(1+2000*exp(-0.6*r)); 
				}				
				r=0.0;
				r=(NYG-jj)*1.0;
				G[1]=0.0;
				G[1]=1/(1+2000*exp(-0.6*r));
				
				grad1[j][i][k]=grad1[j][i][k]*G[0]*G[1];
				grad2[j][i][k]=grad2[j][i][k]*G[0]*G[1];
				grad3[j][i][k]=grad3[j][i][k]*G[0]*G[1];
				
				/*damping x-direction*/
				r=0.0;
				r=(ii)*1.0;
				G[0]=0.0;
				G[0]=1/(1+2000*exp(-0.6*r)); 
				
				r=0.0;
				r=(NXG-ii)*1.0;
				G[1]=0.0;
				G[1]=1/(1+2000*exp(-0.6*r));
				
				grad1[j][i][k]=grad1[j][i][k]*G[0]*G[1];
				grad2[j][i][k]=grad2[j][i][k]*G[0]*G[1];
				grad3[j][i][k]=grad3[j][i][k]*G[0]*G[1];
				
				/*damping z-direction*/
				r=0.0;
				r=(kk)*1.0;
				G[0]=0.0;
				G[0]=1/(1+2000*exp(-0.6*r)); 
				
				r=0.0;
				r=(NZG-kk)*1.0;
				G[1]=0.0;
				G[1]=1/(1+2000*exp(-0.6*r));
				
				grad1[j][i][k]=grad1[j][i][k]*G[0]*G[1];
				grad2[j][i][k]=grad2[j][i][k]*G[0]*G[1];
				grad3[j][i][k]=grad3[j][i][k]*G[0]*G[1];
				/*----------------------------------------------*/
				}
			}
		}
	}



		
}
