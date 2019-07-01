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
 * creates 2.5D model from 3D model, (attention: not in parallel!!); 
 * 2.5D model parameters constant in z-direction
 * S. Butzer 2013
 * ------------------------------------------------------------------------*/

#include "fd.h"

void model2_5(float  ***  rho, float ***  pi, float ***  u, 
float ***  taus, float ***  taup){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern int NXG, NYG, NZG,L;
	/*extern FILE *FP;*/
	
	float ts, tp, muv, Rho, piv;


	int i, j, k;

		for (j=1;j<=NYG;j++){
		      	for (i=1;i<=NXG;i++){
				if(L){
					ts=taus[j][i][100];
					tp=taup[j][i][100];
				}
				muv=u[j][i][100];
				Rho=rho[j][i][100];
				piv=pi[j][i][100];
				for (k=1;k<=NZG;k++){
					if(L){
						taus[j][i][k]=ts;
						taup[j][i][k]=tp;
					}	
					u[j][i][k]=muv;
					rho[j][i][k]=Rho;
					pi[j][i][k]=piv;
				}
			}		
		}
}