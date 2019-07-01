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
*   Calculating source signal at different source positions with different
*   time-shift, centre frequency and amplitude (as specified in SOURCE_FILE).
*   Source signals are written to array signals 
*  ----------------------------------------------------------------------*/

#include "fd.h"


void wavelet(float ** srcpos_loc, int nsrc, int sourceshape, float **signals){


	/* extern variables */
	extern int SOURCE_SHAPE, NT, MYID;
	extern float  DT;
	extern char SIGNAL_FILE[STRING_SIZE];
	extern FILE *FP;

	/*local variables */
	int nts, nt, k;
	float *psource=NULL, tshift, amp=0.0, a, fc, tau, t, ts;


	if (sourceshape==3) psource=rd_sour(&nts,fopen(SIGNAL_FILE,"r"));
	
	/*signals=fmatrix(1,nsrc,1,NT);*/
	
	for (nt=1;nt<=NT;nt++){
			t=(float)nt*DT;
			
			for (k=1;k<=nsrc;k++) {
				tshift=srcpos_loc[4][k];
				fc=srcpos_loc[5][k];
				a=srcpos_loc[6][k];
				ts=1.0/fc;

				switch (sourceshape){
					case 1 : 
						tau=PI*(t-1.5*ts-tshift)/(ts); /* Ricker */
						amp=(((1.0-2.0*tau*tau)*exp(-tau*tau)));
					break;
					case 2 : 
						if ((t<tshift) || (t>(tshift+ts))) amp=0.0;
						else amp=((sin(2.0*PI*(t-tshift)*fc) 
			    				-0.5*sin(4.0*PI*(t-tshift)*fc)));

/*						amp=((sin(2.0*PI*(t+tshift)*fc) 
			    				-0.5*sin(4.0*PI*(t+tshift)*fc)));
*/
					break;
					case 3 : 
						amp=psource[nt];
					break;  /* source wavelet from file SOURCE_FILE */
					case 4 : 
						if ((t<tshift) || (t>(tshift+ts))) amp=0.0;
						else amp=(0.75*PI/ts)*(pow(sin(PI*(t-tshift)/ts),3.0));
						break; /* sinus raised to the power of three */
					case 5: if(t<tshift)amp=0.0;
						else amp=1.0; break;
						
					/*case 6: if(t==tstep) amp=1.0;
						else amp=0.0;
						break;*/
					default : 
						err("Which source-wavelet ? ");
					}
					
					
					signals[k][nt]=amp*a;
		}
	}
	
	fprintf(FP," Message from function wavelet written by PE %d, sourceshape %d\n",MYID, sourceshape);
	fprintf(FP," %d source positions located in subdomain of PE %d \n",nsrc,MYID);
	fprintf(FP," have been assigned with a source signal. \n");
			
		
	if (SOURCE_SHAPE==3) free_vector(psource,1,NT);

}
