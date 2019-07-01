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
 Copying of model parmeters into a testmodel 
 S. Butzer 2013
 *  ----------------------------------------------------------------------*/

#include "fd.h"
void cpmodel(int nx, int ny, int nz, float ***rho, float ***pi, float ***u,float  ***  testrho, float ***  testpi, float ***  testu){

		int j,i,k;
		
		for (j=1;j<=ny;j++){
			for (i=1;i<=nx;i++){
				for (k=1;k<=nz;k++){
				  
				 testrho[j][i][k]=0.0;
				 testrho[j][i][k]=rho[j][i][k];
				 
				 testu[j][i][k]=0.0;
				 testu[j][i][k]=u[j][i][k];
				 
				 testpi[j][i][k]=0.0;
				 testpi[j][i][k]=pi[j][i][k];
				 
				}
			}
		}
		
}