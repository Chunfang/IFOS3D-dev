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
 *   Read Hessian from files  
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void readhess(int nx, int ny, int nz, float ***  hess1, float ***  hess2, float ***hess3, float finv, int iteration){

	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4];
	extern float VP0, VS0, RHO0;
	extern FILE *FP;
	extern char  HESS_FILE[STRING_SIZE];
		
	/* local variables */
	float rhov, vp, vs; /*qp, qs;*/
	int i, j, k, ii, jj, kk;
	FILE *fp_vs, *fp_vp, *fp_rho;
	char filename[STRING_SIZE];

	/* choose data format: ASCII: format=2, BINARY: format=3*/
	const int format=3; 
	

		

	   fprintf(FP,"\n...reading hess information from hess-files...\n");

	   /*sprintf(filename,"hess/hess.vp");*/
	   sprintf(filename,"%s.vp_%4.2fHz_it%d",HESS_FILE,finv,iteration);
	   fp_vp=fopen(filename,"r");
	   if (fp_vp==NULL) err(" Could not open hess_vp ! ");

	   sprintf(filename,"%s.vs_%4.2fHz_it%d",HESS_FILE,finv,iteration);
	   fp_vs=fopen(filename,"r");
	   if (fp_vs==NULL) err(" Could not open  hess_vs! ");

	   sprintf(filename,"%s.rho_%4.2fHz_it%d",HESS_FILE,finv,iteration);
	   fp_rho=fopen(filename,"r");
	   if (fp_rho==NULL) err(" Could not open hess_rho ! ");
	   

	
	/* loop over global grid */
		for (k=1;k<=NZG;k++){
			for (i=1;i<=NXG;i++){
				for (j=1;j<=NYG;j++){

				vp=readdsk(fp_vp, format);
				vs=readdsk(fp_vs, format);
				rhov=readdsk(fp_rho , format);
	
				/* only the PE which belongs to the current global gridpoint 
				is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY)) && 
				    (POS[3]==((k-1)/NZ))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
					kk=k-POS[3]*NZ;
					
					/*hess1[jj][ii][kk]=vp;
					hess2[jj][ii][kk]=vs;
					hess3[jj][ii][kk]=rhov;*/
					
					
					hess1[jj][ii][kk]=vp*pow(VP0,2);
					hess2[jj][ii][kk]=vs*pow(VS0,2);
					hess3[jj][ii][kk]=rhov*pow(RHO0,2);
					
			
				}
				}
			}
		}

	fclose(fp_vp);
	fclose(fp_vs);
	fclose(fp_rho);	
	

}




