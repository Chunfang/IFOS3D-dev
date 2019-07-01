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

/*------------------------------------------------------------
 * steplength calculation using a parabola method:
 * S. Butzer 2013 
 -------------------------------------------------------------*/

#include "fd.h"
void steplength(float *L2, float *step, int iteration, int it_group){

		extern FILE *FP;
		extern FILE *FI;
		extern int MYID, LBFGS;
		extern float TESTSTEP;
		float a=0.0, b=0.0, c=0.0;
		float teststep=0.0;
  
		teststep=TESTSTEP;		

		if(MYID==0)fprintf(FP,"Start estimation of optimal steplength\n");
		
		/*calculation of the steplength parabola*/
		
		a=((step[0]-step[1])*(L2[0]-L2[2])-(step[0]-step[2])*(L2[0]-L2[1]))/(step[0]*step[2]*(step[0]-step[2])+step[0]*step[1]*(step[1]-step[0])+step[1]*step[2]*(step[2]-step[1]));
		b=(L2[0]-L2[1])/(step[0]-step[1]);
		b=b-a*(step[0]+step[1]);
		c=L2[0]-a*step[0]*step[0]-b*step[0];
		
		step[3]=0.0;
		if(a>0){
			step[3]=-b/(2*a);
		}
		
		if(a<0){
			if(L2[0]<L2[1]&&L2[0]<L2[2]) step[3]=0.1*step[1];
			if(L2[1]<L2[0]&&L2[1]<L2[2]) step[3]=step[1];
			if(L2[2]<L2[0]&&L2[2]<L2[1]) step[3]=step[2];
			fprintf(FP," a<0, found maxima\n");
		}
		if(MYID==0)fprintf(FI,"steplength =%e", step[3]);
		
		if(step[3]>2.5*step[4]){ 
			step[3]=2.5*step[4];
		if(MYID==0)fprintf(FP," steplength larger stepmax, set to stepmax\n");
		}
		
		if(LBFGS && it_group>1 && step[3]>1.0){
			step[3]=1.0;
			if(MYID==0)fprintf(FP,"steplength larger 1.0, set to 1.0");
		}
				
		if(step[3]<=0){
			 step[3]=0.1*step[1];
			 if(MYID==0)fprintf(FP," steplength smaller 0,  set to step[3]\n");
		}
		
		if(step[3]<0.1*step[1]){ step[3]=0.1*step[1];
			if(MYID==0)fprintf(FP," Warning: optimal steplength smaller 0.1*step[1]\n");
		}

		L2[3]=a*step[3]*step[3]+b*step[3]+c;
		
		if(MYID==0)fprintf(FP,"\n Steplengthparabel: L2=%.2e*step^2 + %.2e*step + %.2e, Minimum L2(%.2e)=%.2e\n", a,b,c,step[3],L2[3]);
		if(MYID==0)fprintf(FI,"\nSteplengthparabel: L2=%e*step^2+%e*step+%e, Minimum L2(%e)=%e\n", a,b,c,step[3],L2[3]);
		
		/*estimation of new test-steplength*/
		step[4]=step[3]/2;
		
		if(step[4]<0.2*teststep) step[4]=0.2*teststep;
		if(step[4]>teststep) step[4]=teststep;
		if(LBFGS && it_group>1) step[4]=0.5;
		if(MYID==0)fprintf(FI,"new test-steplength4: %e\n", step[4]);

}
