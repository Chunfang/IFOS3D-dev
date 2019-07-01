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
 * 
 * S. Butzer 2015
 * ------------------------------------------------------------------------*/

#include "fd.h"

void constant_boundary(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, float ***  pi, float ***  u, float *** rho){


	/*extern int FREE_SURF;*/
	extern int NPROCX, NPROCY, NPROCZ, POS[4];
	extern int FW,NX, NY;

	int i, j, k,dum1,dum2;
	int dist=4;
	dum1=1;
	dum2=NY;
	

	if (POS[1]==0){	
		if (POS[2]==0){
			dum1=2;
			for (i=1;i<=FW+dist;i++){
				for (k=nz1;k<=nz2;k++){
					u[1][i][k]=u[1][FW+1+dist][k];
					pi[1][i][k]= pi[1][FW+1+dist][k];
					rho[1][i][k]= rho[1][FW+1+dist][k];
				}
			}
		}
		if (POS[2]==NPROCY-1){
			dum2=NY-1; 
			for (i=1;i<=FW+dist;i++){
				for (k=nz1;k<=nz2;k++){
					u[NY][i][k]=u[NY][FW+1+dist][k];
					pi[NY][i][k]= pi[NY][FW+1+dist][k];
					rho[NY][i][k]= rho[NY][FW+1+dist][k];
				}
			}
		}
		for (j=dum1;j<=dum2;j++){
			for (i=1;i<=FW+dist;i++){
				for (k=nz1;k<=nz2;k++){
					u[j][i][k]=(3.0*u[j][FW+1+dist][k]+u[j-1][FW+1+dist][k]+u[j+1][FW+1+dist][k]+u[j][FW+1+dist][k-1]+u[j][FW+1+dist][k+1]+u[j][FW+2+dist][k])/8.0;
					pi[j][i][k]= (3.0*pi[j][FW+1+dist][k]+pi[j-1][FW+1+dist][k]+pi[j+1][FW+1+dist][k]+pi[j][FW+1+dist][k-1]+pi[j][FW+1+dist][k+1]+pi[j][FW+2+dist][k])/8.0;
					rho[j][i][k]= (3.0*rho[j][FW+1+dist][k]+rho[j-1][FW+1+dist][k]+rho[j+1][FW+1+dist][k]+rho[j][FW+1+dist][k-1]+rho[j][FW+1+dist][k+1]+rho[j][FW+2+dist][k])/8.0;
				}
			}
		}

	}
	

	if(POS[1]==NPROCX-1){
	  	if (POS[2]==0){
			dum1=2;
			for (i=nx2+1-dist;i<=nx2+FW;i++){
				for (k=nz1;k<=nz2;k++){
					u[1][i][k]=u[1][nx2-dist][k];
					pi[1][i][k]= pi[1][nx2-dist][k];
					rho[1][i][k]= rho[1][nx2-dist][k];
				}
			}
		}
		if (POS[2]==NPROCY-1){
			dum2=NY-1; 
			for (i=nx2+1-dist;i<=nx2+FW;i++){
				for (k=nz1;k<=nz2;k++){
					u[NY][i][k]=u[NY][nx2-dist][k];
					pi[NY][i][k]= pi[NY][nx2-dist][k];
					rho[NY][i][k]= rho[NY][nx2-dist][k];
				}
			}
		}
		for (j=dum1;j<=dum2;j++){
			for (i=nx2+1-dist;i<=nx2+FW;i++){
				for (k=nz1;k<=nz2;k++){
					u[j][i][k]=(3.0*u[j][nx2-dist][k]+u[j-1][nx2-dist][k]+u[j+1][nx2-dist][k]+u[j][nx2-dist][k-1]+u[j][nx2-dist][k+1]+u[j][nx2-1-dist][k])/8.0;
					pi[j][i][k]=(3.0*pi[j][nx2-dist][k]+pi[j-1][nx2-dist][k]+pi[j+1][nx2-dist][k]+pi[j][nx2-dist][k-1]+pi[j][nx2-dist][k+1]+pi[j][nx2-1-dist][k])/8.0;
					rho[j][i][k]=(3.0*rho[j][nx2-dist][k]+rho[j-1][nx2-dist][k]+rho[j+1][nx2-dist][k]+rho[j][nx2-dist][k-1]+rho[j][nx2-dist][k+1]+rho[j][nx2-1-dist][k])/8.0;
				}
			}
		}
	}

	if(POS[3]==0){
		if (POS[2]==0){
			dum1=2;
			for (i=1;i<=NX;i++){
				for (k=1;k<=FW+dist;k++){
					u[1][i][k]=u[1][i][FW+1+dist];
					pi[1][i][k]= pi[1][i][FW+1+dist];
					rho[1][i][k]= rho[1][i][FW+1+dist];
				}
			}
		}
		if (POS[2]==NPROCY-1){
			dum2=NY-1; 
			for (i=1;i<=NX;i++){
				for (k=1;k<=FW+dist;k++){
					u[NY][i][k]=u[NY][i][FW+1+dist];
					pi[NY][i][k]= pi[NY][i][FW+1+dist];	
					rho[NY][i][k]= rho[NY][i][FW+1+dist];
				}
			}
		}
		if(POS[1]==0){
			for (j=dum1;j<=dum2;j++){
				for (i=1;i<=FW+dist;i++){
					for (k=1;k<=FW+dist;k++){
						u[j][i][k]=u[j][i][FW+1+dist];
						pi[j][i][k]=pi[j][i][FW+1+dist];
						rho[j][i][k]=rho[j][i][FW+1+dist];
					}
				}
			}
		}
		if(POS[1]==NPROCX-1){
			for (j=dum1;j<=dum2;j++){
				for (i=nx2+1-dist;i<=nx2+FW;i++){
					for (k=1;k<=FW+dist;k++){
						u[j][i][k]=u[j][i][FW+1+dist];
						pi[j][i][k]=pi[j][i][FW+1+dist];
						rho[j][i][k]=rho[j][i][FW+1+dist];
					}
				}
			}
		}
		
		
		
		for (j=dum1;j<=dum2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=1;k<=FW+dist;k++){
					u[j][i][k]=(3.0*u[j][i][FW+1+dist]+u[j-1][i][FW+1+dist]+u[j+1][i][FW+1+dist]+u[j][i-1][FW+1+dist]+u[j][i+1][FW+1+dist]+u[j][i][FW+2+dist])/8.0;
					pi[j][i][k]=(3.0*pi[j][i][FW+1+dist]+pi[j-1][i][FW+1+dist]+pi[j+1][i][FW+1+dist]+pi[j][i-1][FW+1+dist]+pi[j][i+1][FW+1+dist]+pi[j][i][FW+2+dist])/8.0;
					rho[j][i][k]=(3.0*rho[j][i][FW+1+dist]+rho[j-1][i][FW+1+dist]+rho[j+1][i][FW+1+dist]+rho[j][i-1][FW+1+dist]+rho[j][i+1][FW+1+dist]+rho[j][i][FW+2+dist])/8.0;
				}
			}
		}
	
	
	}


	if(POS[3]==NPROCZ-1){	
		if (POS[2]==0){
			dum1=2;
			for (i=1;i<=NX;i++){
				for (k=nz2+1-dist;k<=nz2+FW;k++){
					u[1][i][k]=u[1][i][nz2-dist];
					pi[1][i][k]= pi[1][i][nz2-dist];
					rho[1][i][k]= rho[1][i][nz2-dist];
				}
			}
		}
		if (POS[2]==NPROCY-1){
			dum2=NY-1; 
			for (i=1;i<=NX;i++){
				for (k=nz2+1-dist;k<=nz2+FW;k++){
					u[NY][i][k]=u[NY][i][nz2-dist];
					pi[NY][i][k]= pi[NY][i][nz2-dist];
					rho[NY][i][k]= rho[NY][i][nz2-dist];
				}
			}
		}
		if(POS[1]==0){
			for (j=dum1;j<=dum2;j++){
				for (i=1;i<=FW+dist;i++){
					for (k=nz2+1-dist;k<=nz2+FW;k++){
						u[j][i][k]=u[j][i][nz2-dist];
						pi[j][i][k]=pi[j][i][nz2-dist];
						rho[j][i][k]=rho[j][i][nz2-dist];
					}
				}
			}
		}
		if(POS[1]==NPROCX-1){
			for (j=dum1;j<=dum2;j++){
				for (i=nx2+1-dist;i<=nx2+FW;i++){
					for (k=nz2+1-dist;k<=nz2+FW;k++){
						u[j][i][k]=u[j][i][nz2-dist];
						pi[j][i][k]=pi[j][i][nz2-dist];
						rho[j][i][k]=rho[j][i][nz2-dist];
					}
				}
			}
		}
		for (j=dum1;j<=dum2;j++){
			for (i=nx1;i<=nx2;i++){
				for (k=nz2+1-dist;k<=nz2+FW;k++){
					u[j][i][k]=(3.0*u[j][i][nz2-dist]+u[j-1][i][nz2-dist]+u[j+1][i][nz2-dist]+u[j][i-1][nz2-dist]+u[j][i+1][nz2-dist]+u[j][i][nz2-1-dist])/8.0;
					pi[j][i][k]=(3.0*pi[j][i][nz2-dist]+pi[j-1][i][nz2-dist]+pi[j+1][i][nz2-dist]+pi[j][i-1][nz2-dist]+pi[j][i+1][nz2-dist]+pi[j][i][nz2-dist-1])/8.0;
					rho[j][i][k]=(3.0*rho[j][i][nz2-dist]+rho[j-1][i][nz2-dist]+rho[j+1][i][nz2-dist]+rho[j][i-1][nz2-dist]+rho[j][i+1][nz2-dist]+rho[j][i][nz2-dist-1])/8.0;
				}
			}
		}

	}


}
