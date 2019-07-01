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
 * Initialise wavefield (frequency domain) with zero
 -------------------------------------------------------------------------*/

#include "fd.h"

void zero_invers(int NX, int NY, int NZ, float **** Ffvx, float **** Ffvy, float **** Ffvz, float **** Ffivx, float **** Ffivy, float **** Ffivz, float **** Fbvx, float **** Fbvy, float **** Fbvz, float **** Fbivx, float **** Fbivy, float **** Fbivz,int nfmax, int ntr_hess){

	extern int FDORDER, ABS_TYPE, POS[4];
	int nx1, ny1, nz1, nx2, ny2, nz2, a,b,l, i, j, k,m;
	
	
	l=1;
	if(ABS_TYPE==1 && FDORDER==2){l=2;}
	    
	    if(POS[2]==0){
	      a=0;
	      b=1;}
	    else{
	      a=1;
	      b=l;}
	      

	ny1=a-b*FDORDER/2;
	ny2=NY+l*FDORDER/2;
	nx1=1-l*FDORDER/2;
	nx2=NX+l*FDORDER/2;
	nz1=1-l*FDORDER/2;
	nz2=NZ+l*FDORDER/2;

		
	for(m=1;m<=nfmax*(ntr_hess+1);m++){
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){
				  Fbvx[m][j][i][k]=0.0;
				  Fbvy[m][j][i][k]=0.0;
				  Fbvz[m][j][i][k]=0.0;
				  Fbivx[m][j][i][k]=0.0;
				  Fbivy[m][j][i][k]=0.0;
				  Fbivz[m][j][i][k]=0.0;
				}
			}
		      }
	    }
	    
	    	for(m=1;m<=nfmax;m++){
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){
				  Ffvx[m][j][i][k]=0.0;
				  Ffvy[m][j][i][k]=0.0;
				  Ffvz[m][j][i][k]=0.0;
				  Ffivx[m][j][i][k]=0.0;
				  Ffivy[m][j][i][k]=0.0;
				  Ffivz[m][j][i][k]=0.0;
				}
			}
		      }
	    }

}