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
 * calculation of diagonal Hessian approximation in frequency domain:
 * spatial derivatives are calculated by 4th order finite differences
 * S. Butzer 2013
--------------------------------------------------------------------------- */

#include "fd.h"

void hess_F(int nx,int ny,int nz,float **** fvx,float **** fvy,float **** fvz,float **** fivx,float **** fivy,float **** fivz,float **** bvx,float **** bvy,float **** bvz, float **** bivx,float **** bivy,float **** bivz, float ***hess1, float ***hess2,float ***hess3,int nt, float  ***  rho, float ***  pi, float ***  u, float * finv, int nf, int ntr_hess){

	extern float DX, DY, DZ;
	extern int  FDCOEFF;
	/*extern char  MFILE[STRING_SIZE];*/
		
	float fvxx=0.0,fvxy=0.0,fvxz=0.0,fvyx=0.0,fvyy=0.0,fvyz=0.0,fvzx=0.0,fvzy=0.0,fvzz=0.0;
	float bvxx=0.0,bvxy=0.0,bvxz=0.0,bvyx=0.0,bvyy=0.0,bvyz=0.0,bvzx=0.0,bvzy=0.0,bvzz=0.0;
	float fivxx=0.0,fivxy=0.0,fivxz=0.0,fivyx=0.0,fivyy=0.0,fivyz=0.0,fivzx=0.0,fivzy=0.0,fivzz=0.0;
	float bivxx=0.0,bivxy=0.0,bivxz=0.0,bivyx=0.0,bivyy=0.0,bivyz=0.0,bivzx=0.0,bivzy=0.0,bivzz=0.0;
	float relam, remu, rerho,imlam, immu, imrho,revs,imvs, rerho1, imrho1;
	/*float hessmu=0.0, hesslam=0.0, hessrho=0.0;*/
	float b1,b2,fdummy;
	/*float vp0=6200.0, vs0=3600.0, rho0=2800.0;*/
	
	int i,j,k,l,m,n;
	/*char gradfile1[STRING_SIZE],gradfile2[STRING_SIZE],gradfile3[STRING_SIZE],gradfile4[STRING_SIZE],gradfile5[STRING_SIZE],gradfile6[STRING_SIZE],gradfile7[STRING_SIZE];
	FILE *fpmod1, *fpmod2, *fpmod3,*fpmod4, *fpmod5, *fpmod6,*fpmod7;
	
	sprintf(gradfile1,"%s.grad1.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile2,"%s.grad2.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile3,"%s.grad3.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile4,"%s.grad4.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile5,"%s.grad5.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile6,"%s.grad6.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	sprintf(gradfile7,"%s.grad7.%i%i%i",MFILE,POS[1],POS[2],POS[3]);
	
	fpmod1=fopen(gradfile1,"w");
	fpmod2=fopen(gradfile2,"w");
	fpmod3=fopen(gradfile3,"w");
	fpmod4=fopen(gradfile4,"w");
	fpmod5=fopen(gradfile5,"w");
	fpmod6=fopen(gradfile6,"w");
	fpmod7=fopen(gradfile7,"w");*/
	
	
	/*t=NT-nt+1;*/
	b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients (4th order)*/
		if(FDCOEFF==2){
		b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/	
	for(n=1;n<=ntr_hess;n++)
	for(l=1;l<=nf;l++){
	  m=(n)*nf+l;
	  fdummy=0.0;
	  fdummy=finv[l-1]*M_PI*2;
	  fdummy=1/fdummy;
		for (j=1;j<=ny;j++){
			for (i=1;i<=nx;i++){
				for (k=1;k<=nz;k++){

				/* spatial derivatives of the components of the velocities
						    are computed */
								    
				fvxx = (b1*(fvx[l][j][i][k]-fvx[l][j][i-1][k])+b2*(fvx[l][j][i+1][k]-fvx[l][j][i-2][k]))/DX;
				fvxy = (b1*(fvx[l][j+1][i][k]-fvx[l][j][i][k])+b2*(fvx[l][j+2][i][k]-fvx[l][j-1][i][k]))/DY;
				fvxz = (b1*(fvx[l][j][i][k+1]-fvx[l][j][i][k])+b2*(fvx[l][j][i][k+2]-fvx[l][j][i][k-1]))/DZ;		    
				fvyx = (b1*(fvy[l][j][i+1][k]-fvy[l][j][i][k])+b2*(fvy[l][j][i+2][k]-fvy[l][j][i-1][k]))/DX;
                                fvyy = (b1*(fvy[l][j][i][k]-fvy[l][j-1][i][k])+b2*(fvy[l][j+1][i][k]-fvy[l][j-2][i][k]))/DY;
				fvyz = (b1*(fvy[l][j][i][k+1]-fvy[l][j][i][k])+b2*(fvy[l][j][i][k+2]-fvy[l][j][i][k-1]))/DZ;
			        fvzx = (b1*(fvz[l][j][i+1][k]-fvz[l][j][i][k])+b2*(fvz[l][j][i+2][k]-fvz[l][j][i-1][k]))/DX;
				fvzy = (b1*(fvz[l][j+1][i][k]-fvz[l][j][i][k])+b2*(fvz[l][j+2][i][k]-fvz[l][j-1][i][k]))/DY;
				fvzz = (b1*(fvz[l][j][i][k]-fvz[l][j][i][k-1])+b2*(fvz[l][j][i][k+1]-fvz[l][j][i][k-2]))/DZ;

				bvxx = (b1*(bvx[m][j][i][k]-bvx[m][j][i-1][k])+b2*(bvx[m][j][i+1][k]-bvx[m][j][i-2][k]))/DX;
				bvxy = (b1*(bvx[m][j+1][i][k]-bvx[m][j][i][k])+b2*(bvx[m][j+2][i][k]-bvx[m][j-1][i][k]))/DY;		    
         			bvxz = (b1*(bvx[m][j][i][k+1]-bvx[m][j][i][k])+b2*(bvx[m][j][i][k+2]-bvx[m][j][i][k-1]))/DZ;		    
				bvyx = (b1*(bvy[m][j][i+1][k]-bvy[m][j][i][k])+b2*(bvy[m][j][i+2][k]-bvy[m][j][i-1][k]))/DX;
                                bvyy = (b1*(bvy[m][j][i][k]-bvy[m][j-1][i][k])+b2*(bvy[m][j+1][i][k]-bvy[m][j-2][i][k]))/DY;
				bvyz = (b1*(bvy[m][j][i][k+1]-bvy[m][j][i][k])+b2*(bvy[m][j][i][k+2]-bvy[m][j][i][k-1]))/DZ;
			        bvzx = (b1*(bvz[m][j][i+1][k]-bvz[m][j][i][k])+b2*(bvz[m][j][i+2][k]-bvz[m][j][i-1][k]))/DX;
				bvzy = (b1*(bvz[m][j+1][i][k]-bvz[m][j][i][k])+b2*(bvz[m][j+2][i][k]-bvz[m][j-1][i][k]))/DY;
				bvzz = (b1*(bvz[m][j][i][k]-bvz[m][j][i][k-1])+b2*(bvz[m][j][i][k+1]-bvz[m][j][i][k-2]))/DZ;
				
				fivxx = (b1*(fivx[l][j][i][k]-fivx[l][j][i-1][k])+b2*(fivx[l][j][i+1][k]-fivx[l][j][i-2][k]))/DX;
				fivxy = (b1*(fivx[l][j+1][i][k]-fivx[l][j][i][k])+b2*(fivx[l][j+2][i][k]-fivx[l][j-1][i][k]))/DY;		    
         			fivxz = (b1*(fivx[l][j][i][k+1]-fivx[l][j][i][k])+b2*(fivx[l][j][i][k+2]-fivx[l][j][i][k-1]))/DZ;		    
				fivyx = (b1*(fivy[l][j][i+1][k]-fivy[l][j][i][k])+b2*(fivy[l][j][i+2][k]-fivy[l][j][i-1][k]))/DX;
                                fivyy = (b1*(fivy[l][j][i][k]-fivy[l][j-1][i][k])+b2*(fivy[l][j+1][i][k]-fivy[l][j-2][i][k]))/DY;
				fivyz = (b1*(fivy[l][j][i][k+1]-fivy[l][j][i][k])+b2*(fivy[l][j][i][k+2]-fivy[l][j][i][k-1]))/DZ;
			        fivzx = (b1*(fivz[l][j][i+1][k]-fivz[l][j][i][k])+b2*(fivz[l][j][i+2][k]-fivz[l][j][i-1][k]))/DX;
				fivzy = (b1*(fivz[l][j+1][i][k]-fivz[l][j][i][k])+b2*(fivz[l][j+2][i][k]-fivz[l][j-1][i][k]))/DY;
				fivzz = (b1*(fivz[l][j][i][k]-fivz[l][j][i][k-1])+b2*(fivz[l][j][i][k+1]-fivz[l][j][i][k-2]))/DZ;

				bivxx = (b1*(bivx[m][j][i][k]-bivx[m][j][i-1][k])+b2*(bivx[m][j][i+1][k]-bivx[m][j][i-2][k]))/DX;
				bivxy = (b1*(bivx[m][j+1][i][k]-bivx[m][j][i][k])+b2*(bivx[m][j+2][i][k]-bivx[m][j-1][i][k]))/DY;		    
         			bivxz = (b1*(bivx[m][j][i][k+1]-bivx[m][j][i][k])+b2*(bivx[m][j][i][k+2]-bivx[m][j][i][k-1]))/DZ;		    
				bivyx = (b1*(bivy[m][j][i+1][k]-bivy[m][j][i][k])+b2*(bivy[m][j][i+2][k]-bivy[m][j][i-1][k]))/DX;
                                bivyy = (b1*(bivy[m][j][i][k]-bivy[m][j-1][i][k])+b2*(bivy[m][j+1][i][k]-bivy[m][j-2][i][k]))/DY;
				bivyz = (b1*(bivy[m][j][i][k+1]-bivy[m][j][i][k])+b2*(bivy[m][j][i][k+2]-bivy[m][j][i][k-1]))/DZ;
			        bivzx = (b1*(bivz[m][j][i+1][k]-bivz[m][j][i][k])+b2*(bivz[m][j][i+2][k]-bivz[m][j][i-1][k]))/DX;
				bivzy = (b1*(bivz[m][j+1][i][k]-bivz[m][j][i][k])+b2*(bivz[m][j+2][i][k]-bivz[m][j-1][i][k]))/DY;
				bivzz = (b1*(bivz[m][j][i][k]-bivz[m][j][i][k-1])+b2*(bivz[m][j][i][k+1]-bivz[m][j][i][k-2]))/DZ;
				
				
			

				relam=0.0; imlam=0.0;   /*relam and imlam correponds to partial derivative wavefields with respect to lambda*/ /*Geändert wg. f als Geschwindigkeit, b als displacement*/
				relam=(fivxx+fivyy+fivzz)*(bvxx+bvyy+bvzz)+(fvxx+fvyy+fvzz)*(bivxx+bivyy+bivzz);
				relam=fdummy*relam;
				imlam=(fivxx+fivyy+fivzz)*(bivxx+bivyy+bivzz)-(fvxx+fvyy+fvzz)*(bvxx+bvyy+bvzz);
				imlam=fdummy*imlam;
				
				remu=0.0; immu=0.0;		
				remu=2*fivxx*bvxx+2*fivyy*bvyy+2*fivzz*bvzz+2*fvxx*bivxx+2*fvyy*bivyy+2*fvzz*bivzz+(fivxy+fivyx)*(bvxy+bvyx)+(fvxy+fvyx)*(bivxy+bivyx)+(fivxz+fivzx)*(bvxz+bvzx)+(fvxz+fvzx)*(bivxz+bivzx)+(fivyz+fivzy)*(bvyz+bvzy)+(fvyz+fvzy)*(bivyz+bivzy);
				remu=fdummy*remu;	
				immu=2*fivxx*bivxx+2*fivyy*bivyy+2*fivzz*bivzz-2*fvxx*bvxx-2*fvyy*bvyy-2*fvzz*bvzz+(fivxy+fivyx)*(bivxy+bivyx)-(fvxy+fvyx)*(bvxy+bvyx)+(fivxz+fivzx)*(bivxz+bivzx)-(fvxz+fvzx)*(bvxz+bvzx)+(fivyz+fivzy)*(bivyz+bivzy)-(fvyz+fvzy)*(bvyz+bvzy);
				immu=fdummy*immu;
				
				rerho=0.0; imrho=0.0;/*rho noch nicht abgeändert!!*/
				rerho=(fvx[l][j][i][k]*bvx[m][j][i][k]+fvy[l][j][i][k]*bvy[m][j][i][k]+fvz[l][j][i][k]*bvz[m][j][i][k])- (fivx[l][j][i][k]*bivx[m][j][i][k]+fivy[l][j][i][k]*bivy[m][j][i][k]+fivz[l][j][i][k]*bivz[m][j][i][k]);
				imrho=(fvx[l][j][i][k]*bivx[m][j][i][k]+fvy[l][j][i][k]*bivy[m][j][i][k]+fvz[l][j][i][k]*bivz[m][j][i][k])+ (fivx[l][j][i][k]*bvx[m][j][i][k]+fivy[l][j][i][k]*bvy[m][j][i][k]+fivz[l][j][i][k]*bvz[m][j][i][k]);
				
				revs=0.0; imvs=0.0;
				revs=-2*relam+remu;
				imvs=-2*imlam+immu;
				
				rerho1=0.0; imrho1=0.0;
				rerho1=rerho+relam*(pi[j][i][k]-2*u[j][i][k])/rho[j][i][k]+remu*u[j][i][k]/rho[j][i][k];
				imrho1=imrho+imlam*(pi[j][i][k]-2*u[j][i][k])/rho[j][i][k]+immu*u[j][i][k]/rho[j][i][k];
				
				/*hesslam=relam*relam+imlam*imlam;
				hessmu=remu*remu+immu*immu;
				hessrho=rerho*rerho+imrho*imrho;*/
				
				
				hess1[j][i][k]+=(relam*relam+imlam*imlam)*rho[j][i][k]*pi[j][i][k]*4;
				hess2[j][i][k]+=(revs*revs+imvs*imvs)*4*rho[j][i][k]*u[j][i][k];
				hess3[j][i][k]+=rerho1*rerho1+imrho1*imrho1;
				
							
				
				/*fwrite(&hess1[j][i][k], sizeof(float), 1,fpmod1);
				fwrite(&hess2[j][i][k], sizeof(float), 1,fpmod2);
				fwrite(&hess3[j][i][k], sizeof(float), 1,fpmod3);
				fwrite(&bivz[l][j][i][k], sizeof(float), 1,fpmod4);
				fwrite(&bvx[l][j][i][k], sizeof(float), 1,fpmod5);
				fwrite(&fivx[l][j][i][k], sizeof(float), 1,fpmod6);
				fwrite(&fvxy, sizeof(float), 1,fpmod7);*/
	
				}
			}
		}
	}
	
		/*Normierung der Hessian, falls LBFGS*/
	/*for (j=1;j<=ny;j++){
		for (i=1;i<=nx;i++){
			for (k=1;k<=nz;k++){
			hess1[j][i][k]=hess1[j][i][k]*pow(vp0,2.0);
			hess2[j][i][k]=hess2[j][i][k]*pow(vs0,2.0);
			hess3[j][i][k]=hess3[j][i][k]*pow(rho0,2.0);

			}
		}
	}*/

	/*fclose(fpmod1);
	fclose(fpmod2);
	fclose(fpmod3);
	fclose(fpmod4);
	fclose(fpmod5);
	fclose(fpmod6);
	fclose(fpmod7);*/
	

}