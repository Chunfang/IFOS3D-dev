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

/*-------------------------------------------------------------------------
 * calculation of a discrete Fourier transfomation on the fly:
 * Fouriercomponents of forward or backpropagated wavefields are summed up for each frequency
 * S. Butezr 2013
 --------------------------------------------------------------------------*/
 

#include "fd.h"
void discfourier(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, int nt, float *** vx, float *** vy, float *** vz, float **** Fvx, float **** Fvy,float **** Fvz, float **** Fvxi, float **** Fvyi,float **** Fvzi,float *finv, int nf, int ntast, int pshot1, int back){
  
	extern float DT,TIME;
	
	int i, j, k, l,m;
	double trig1,trig2;
	float t=0.0;
  
	if(back==0) t=nt*DT;
	if(back==1) t=TIME-nt*DT;
	
	for(l=1;l<=nf;l++){
		m=(pshot1)*nf+l;
		trig1=0.0;
		trig1=cos(2.0*t*finv[l-1]*M_PI)*DT*ntast;
		trig2=0.0;
		trig2=sin(2.0*t*finv[l-1]*M_PI)*DT*ntast;
	  
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz1;k<=nz2;k++){
				  
					  Fvx[m][j][i][k]+=vx[j][i][k]*trig1;
					  Fvy[m][j][i][k]+=vy[j][i][k]*trig1;
					  Fvz[m][j][i][k]+=vz[j][i][k]*trig1;
					  
					  Fvxi[m][j][i][k]+=vx[j][i][k]*trig2;
					  Fvyi[m][j][i][k]+=vy[j][i][k]*trig2;
					  Fvzi[m][j][i][k]+=vz[j][i][k]*trig2;
		  
				}
			}
		}
	}				
}