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
* Filter seismograms in time domain with a Butterworth filter of the libseife
* library. Lowpass or highpass filtering can be applied                                 
* last update 2011, L. Rehor
*  ----------------------------------------------------------------------*/
	

#include "fd.h"
#include "cseife.h"

void filt_seis(float ** data,int ntr, int ns, float finv){
	  
	  int order=4;
	  int method=1;
	  float fc=0.0;
	  extern int VERBOSE;
	 extern FILE *FP;
	        /*
	        data    :       2-dimensional array containing seismograms (
	        finv    :       corner frequency in Hz
	        order   :       order of filter
	        ntr     :       number of traces
	        ns      :       number of samples
	        method  :       definition of filter
	                        1: lowpass filter
	                        2: highpass filter
	        */
	                  
	
	        /* declaration of extern variables */
	        extern float DT;
	       
	        /* declaration of local variables */
	        int itr, j;
	        double *seismogram, T0;
	       
       
	       fc=finv;
	        seismogram=dvector(1,ns);
	       
	        T0=1.0/fc;
		
		
		if (VERBOSE) fprintf(FP,"  ns=%d; DT=%e; T0=%e",ns,DT,T0);
	       
	        for (itr=1;itr<=ntr;itr++){
	                for (j=1;j<=ns;j++){
	                        seismogram[j]=(double)data[itr][j];}
              
	                if (method==1){         /*lowpass filter*/
	                        seife_lpb(seismogram,ns,DT,T0,order);}
	               
	               
                if (method==2){         /*highpass filter*/
	                        seife_hpb(seismogram,ns,DT,T0,order);}
               
                for (j=1;j<=ns;j++){
	                        data[itr][j]=(float)seismogram[j];}
	        }
	       
	        free_dvector(seismogram,1,ns);
}