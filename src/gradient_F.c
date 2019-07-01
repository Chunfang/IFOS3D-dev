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
 * gradient calculation in frequency domain:
 * gradient as multiplication of forward and conjugate backpropagated wavefield 
 * spatial derivatives are calculated by 4th order finite differences
 * S. Butzer 2013
--------------------------------------------------------------------------- */


#include "fd.h"

void gradient_F(int nx,int ny,int nz,float **** fvx,float **** fvy,float **** fvz,float **** fivx,float **** fivy,float **** fivz,float **** bvx,float **** bvy,float **** bvz, float **** bivx,float **** bivy,float **** bivz, float ***grad1, float ***grad2,float ***grad3,int nt, float  ***  rho, float ***  pi, float ***  u, float * finv, int nf, int iteration){

	extern float DX, DY, DZ;
	extern int POS[4], FDCOEFF;
	extern char  MFILE[STRING_SIZE];
		
	float fvxx=0.0,fvxy=0.0,fvxz=0.0,fvyx=0.0,fvyy=0.0,fvyz=0.0,fvzx=0.0,fvzy=0.0,fvzz=0.0;
	float bvxx=0.0,bvxy=0.0,bvxz=0.0,bvyx=0.0,bvyy=0.0,bvyz=0.0,bvzx=0.0,bvzy=0.0,bvzz=0.0;
	float fivxx=0.0,fivxy=0.0,fivxz=0.0,fivyx=0.0,fivyy=0.0,fivyz=0.0,fivzx=0.0,fivzy=0.0,fivzz=0.0;
	float bivxx=0.0,bivxy=0.0,bivxz=0.0,bivyx=0.0,bivyy=0.0,bivyz=0.0,bivzx=0.0,bivzy=0.0,bivzz=0.0;
	float gradlam, gradmu, gradrho;
	float b1,b2,fdummy;
	/*float vp0=6200.0, vs0=3600.0, rho0=2800.0;*/
	int i,j,k,l;
	char gradfile1[STRING_SIZE],gradfile2[STRING_SIZE],gradfile3[STRING_SIZE],gradfile4[STRING_SIZE],gradfile5[STRING_SIZE],gradfile6[STRING_SIZE],gradfile7[STRING_SIZE];
	
	/*FILE *fpmod1, *fpmod2, *fpmod3,*fpmod4, *fpmod5, *fpmod6,*fpmod7;*/
	
		
	sprintf(gradfile1,"%s.grad1.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile2,"%s.grad2.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile3,"%s.grad3.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile4,"%s.grad4.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile5,"%s.grad5.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile6,"%s.grad6.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile7,"%s.grad7.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	
	/*fpmod1=fopen(gradfile1,"w");
	fpmod2=fopen(gradfile2,"w");
	fpmod3=fopen(gradfile3,"w");
	fpmod4=fopen(gradfile4,"w");
	fpmod5=fopen(gradfile5,"w");
	fpmod6=fopen(gradfile6,"w");
	fpmod7=fopen(gradfile7,"w");*/
	

	b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients (4th order)*/
		if(FDCOEFF==2){
		b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/	
	for(l=1;l<=nf;l++){
		fdummy=0.0;
		fdummy=2.0*2.0*finv[l-1]*finv[l-1]*M_PI*M_PI;
		fdummy=1/fdummy;
		for (j=1;j<=ny;j++){
			for (i=1;i<=nx;i++){
				for (k=1;k<=nz;k++){

				/* spatial derivatives of the components of the velocities are computed */
								    
				fvxx = (b1*(fvx[l][j][i][k]-fvx[l][j][i-1][k])+b2*(fvx[l][j][i+1][k]-fvx[l][j][i-2][k]))/DX;
				fvxy = (b1*(fvx[l][j+1][i][k]-fvx[l][j][i][k])+b2*(fvx[l][j+2][i][k]-fvx[l][j-1][i][k]))/DY;
				fvxz = (b1*(fvx[l][j][i][k+1]-fvx[l][j][i][k])+b2*(fvx[l][j][i][k+2]-fvx[l][j][i][k-1]))/DZ;		    
				fvyx = (b1*(fvy[l][j][i+1][k]-fvy[l][j][i][k])+b2*(fvy[l][j][i+2][k]-fvy[l][j][i-1][k]))/DX;
                                fvyy = (b1*(fvy[l][j][i][k]-fvy[l][j-1][i][k])+b2*(fvy[l][j+1][i][k]-fvy[l][j-2][i][k]))/DY;
				fvyz = (b1*(fvy[l][j][i][k+1]-fvy[l][j][i][k])+b2*(fvy[l][j][i][k+2]-fvy[l][j][i][k-1]))/DZ;
			        fvzx = (b1*(fvz[l][j][i+1][k]-fvz[l][j][i][k])+b2*(fvz[l][j][i+2][k]-fvz[l][j][i-1][k]))/DX;
				fvzy = (b1*(fvz[l][j+1][i][k]-fvz[l][j][i][k])+b2*(fvz[l][j+2][i][k]-fvz[l][j-1][i][k]))/DY;
				fvzz = (b1*(fvz[l][j][i][k]-fvz[l][j][i][k-1])+b2*(fvz[l][j][i][k+1]-fvz[l][j][i][k-2]))/DZ;

				bvxx = (b1*(bvx[l][j][i][k]-bvx[l][j][i-1][k])+b2*(bvx[l][j][i+1][k]-bvx[l][j][i-2][k]))/DX;
				bvxy = (b1*(bvx[l][j+1][i][k]-bvx[l][j][i][k])+b2*(bvx[l][j+2][i][k]-bvx[l][j-1][i][k]))/DY;		    
         			bvxz = (b1*(bvx[l][j][i][k+1]-bvx[l][j][i][k])+b2*(bvx[l][j][i][k+2]-bvx[l][j][i][k-1]))/DZ;		    
				bvyx = (b1*(bvy[l][j][i+1][k]-bvy[l][j][i][k])+b2*(bvy[l][j][i+2][k]-bvy[l][j][i-1][k]))/DX;
                                bvyy = (b1*(bvy[l][j][i][k]-bvy[l][j-1][i][k])+b2*(bvy[l][j+1][i][k]-bvy[l][j-2][i][k]))/DY;
				bvyz = (b1*(bvy[l][j][i][k+1]-bvy[l][j][i][k])+b2*(bvy[l][j][i][k+2]-bvy[l][j][i][k-1]))/DZ;
			        bvzx = (b1*(bvz[l][j][i+1][k]-bvz[l][j][i][k])+b2*(bvz[l][j][i+2][k]-bvz[l][j][i-1][k]))/DX;
				bvzy = (b1*(bvz[l][j+1][i][k]-bvz[l][j][i][k])+b2*(bvz[l][j+2][i][k]-bvz[l][j-1][i][k]))/DY;
				bvzz = (b1*(bvz[l][j][i][k]-bvz[l][j][i][k-1])+b2*(bvz[l][j][i][k+1]-bvz[l][j][i][k-2]))/DZ;
				
				fivxx = (b1*(fivx[l][j][i][k]-fivx[l][j][i-1][k])+b2*(fivx[l][j][i+1][k]-fivx[l][j][i-2][k]))/DX;
				fivxy = (b1*(fivx[l][j+1][i][k]-fivx[l][j][i][k])+b2*(fivx[l][j+2][i][k]-fivx[l][j-1][i][k]))/DY;		    
         			fivxz = (b1*(fivx[l][j][i][k+1]-fivx[l][j][i][k])+b2*(fivx[l][j][i][k+2]-fivx[l][j][i][k-1]))/DZ;		    
				fivyx = (b1*(fivy[l][j][i+1][k]-fivy[l][j][i][k])+b2*(fivy[l][j][i+2][k]-fivy[l][j][i-1][k]))/DX;
                                fivyy = (b1*(fivy[l][j][i][k]-fivy[l][j-1][i][k])+b2*(fivy[l][j+1][i][k]-fivy[l][j-2][i][k]))/DY;
				fivyz = (b1*(fivy[l][j][i][k+1]-fivy[l][j][i][k])+b2*(fivy[l][j][i][k+2]-fivy[l][j][i][k-1]))/DZ;
			        fivzx = (b1*(fivz[l][j][i+1][k]-fivz[l][j][i][k])+b2*(fivz[l][j][i+2][k]-fivz[l][j][i-1][k]))/DX;
				fivzy = (b1*(fivz[l][j+1][i][k]-fivz[l][j][i][k])+b2*(fivz[l][j+2][i][k]-fivz[l][j-1][i][k]))/DY;
				fivzz = (b1*(fivz[l][j][i][k]-fivz[l][j][i][k-1])+b2*(fivz[l][j][i][k+1]-fivz[l][j][i][k-2]))/DZ;

				bivxx = (b1*(bivx[l][j][i][k]-bivx[l][j][i-1][k])+b2*(bivx[l][j][i+1][k]-bivx[l][j][i-2][k]))/DX;
				bivxy = (b1*(bivx[l][j+1][i][k]-bivx[l][j][i][k])+b2*(bivx[l][j+2][i][k]-bivx[l][j-1][i][k]))/DY;		    
         			bivxz = (b1*(bivx[l][j][i][k+1]-bivx[l][j][i][k])+b2*(bivx[l][j][i][k+2]-bivx[l][j][i][k-1]))/DZ;		    
				bivyx = (b1*(bivy[l][j][i+1][k]-bivy[l][j][i][k])+b2*(bivy[l][j][i+2][k]-bivy[l][j][i-1][k]))/DX;
                                bivyy = (b1*(bivy[l][j][i][k]-bivy[l][j-1][i][k])+b2*(bivy[l][j+1][i][k]-bivy[l][j-2][i][k]))/DY;
				bivyz = (b1*(bivy[l][j][i][k+1]-bivy[l][j][i][k])+b2*(bivy[l][j][i][k+2]-bivy[l][j][i][k-1]))/DZ;
			        bivzx = (b1*(bivz[l][j][i+1][k]-bivz[l][j][i][k])+b2*(bivz[l][j][i+2][k]-bivz[l][j][i-1][k]))/DX;
				bivzy = (b1*(bivz[l][j+1][i][k]-bivz[l][j][i][k])+b2*(bivz[l][j+2][i][k]-bivz[l][j-1][i][k]))/DY;
				bivzz = (b1*(bivz[l][j][i][k]-bivz[l][j][i][k-1])+b2*(bivz[l][j][i][k+1]-bivz[l][j][i][k-2]))/DZ;
				
				
			

				gradlam=0.0;
				gradlam=(fvxx+fvyy+fvzz)*(bvxx+bvyy+bvzz)+(fivxx+fivyy+fivzz)*(bivxx+bivyy+bivzz);
				gradlam=-gradlam*fdummy;
				
				gradmu=0.0;
				gradmu= 2*fvxx*bvxx+2*fvyy*bvyy+2*fvzz*bvzz+2*fivxx*bivxx+2*fivyy*bivyy+2*fivzz*bivzz+(fvxy+fvyx)*(bvxy+bvyx)+(fivxy+fivyx)*(bivxy+bivyx)+(fvxz+fvzx)*(bvxz+bvzx)+(fivxz+fivzx)*(bivxz+bivzx)+(fvyz+fvzy)*(bvyz+bvzy)+(fivyz+fivzy)*(bivyz+bivzy);
				gradmu=-gradmu*fdummy;
				
				gradrho=0.0;
				gradrho=(fvx[l][j][i][k]*bvx[l][j][i][k]+fvy[l][j][i][k]*bvy[l][j][i][k]+fvz[l][j][i][k]*bvz[l][j][i][k])+ (fivx[l][j][i][k]*bivx[l][j][i][k]+fivy[l][j][i][k]*bivy[l][j][i][k]+fivz[l][j][i][k]*bivz[l][j][i][k]);
				gradrho=gradrho;
				
				/*parametrisation vp, vs, rho*/
				grad1[j][i][k]+=sqrt(rho[j][i][k]*pi[j][i][k])*2*gradlam; /*gradient vp*/
				
				grad2[j][i][k]+=-4*sqrt(rho[j][i][k]*u[j][i][k])*gradlam+2*sqrt(rho[j][i][k]*u[j][i][k])*gradmu; /*gradient vs*/
		
				grad3[j][i][k]+=gradrho+u[j][i][k]/rho[j][i][k]*gradmu+(pi[j][i][k]-2*u[j][i][k])/rho[j][i][k]*gradlam; /*gradient rho*/
				
							
				/*parametrisation lambda, mu, rho*/
				/*grad1[j][i][k]+=gradlam;  
				grad2[j][i][k]+=gradmu;
				grad3[j][i][k]+=gradrho;*/
				
				/*fwrite(&fvz[l][j][i][k], sizeof(float), 1,fpmod1);
				fwrite(&fivz[l][j][i][k], sizeof(float), 1,fpmod2);
				fwrite(&bvz[l][j][i][k], sizeof(float), 1,fpmod3);
				fwrite(&bivz[l][j][i][k], sizeof(float), 1,fpmod4);
				fwrite(&bvx[l][j][i][k], sizeof(float), 1,fpmod5);
				fwrite(&fivx[l][j][i][k], sizeof(float), 1,fpmod6);
				fwrite(&fvxy, sizeof(float), 1,fpmod7);*/
	
				}
			}
		}
	}
	


	/*fclose(fpmod1);
	fclose(fpmod2);
	fclose(fpmod3);
	fclose(fpmod4);
	fclose(fpmod5);
	fclose(fpmod6);
	fclose(fpmod7);*/
	

}