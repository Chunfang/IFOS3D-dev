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

/*---------------------------------------------------------------------------
 * calculation of data residuals (L2) and L2 norm in reverse time direction 
 --------------------------------------------------------------------------*/


#include "fd.h"

void residual(float **sectiondata, float **sectiondataf, float **section, float **sectiondiff, int ntr, int ns, float *L2, float *L2f ){
    
	extern float DT;
	float vdiff;
	int i,j;
	int rt=0;
	extern float F_INV;
	float f_inv=0.0,t=0.0;
	float sectiondifff=0.0,sectiondifffi=0.0;

	f_inv=F_INV;
	
	if(rt==0){
		for(i=1;i<=ntr;i++){
			vdiff=0.0;
		    
			for(j=1;j<=ns;j++){
				sectiondiff[i][ns-j+1]=0.0;
			      
				vdiff+=section[i][j]-sectiondata[i][j];
				sectiondiff[i][ns-j+1]=DT*vdiff;     /* as displacement */
				t=0.0;
				t=DT*(ns-j+1);
				sectiondifff+=sectiondiff[i][ns-j+1]*cos(2.0*t*f_inv*M_PI)*DT;
				sectiondifffi+=sectiondiff[i][ns-j+1]*sin(2.0*t*f_inv*M_PI)*DT;

				*L2+=fabs(sectiondiff[i][ns-j+1]*sectiondiff[i][ns-j+1]);
		    
		    
			}
			
			*L2f+=sectiondifff*sectiondifff+sectiondifffi*sectiondifffi;
		}
	}
    /*printf("L2=%e", L2);*/
	else{
		for(i=1;i<=ntr;i++){
			vdiff=0.0;
		    
			for(j=1;j<=ns;j++){
			/*	sectiondiff[i][ns-j+1]=0.0;
			      
				vdiff+=sectiondata[i][j];
				sectiondiff[i][ns-j+1]=DT*vdiff; */    /* as displacement */
				
				sectiondiff[i][j]=0.0;
				vdiff+=section[i][j];
			        sectiondiff[i][j]=DT*vdiff;
				/*vdiff+=sectiondata[i][j];
				sectiondiff[i][j]=DT*vdiff;  */ /* as displacement */
				
				
				
			}
		} 
	}


}