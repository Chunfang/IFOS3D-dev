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
 *   smooth model parameters (complete model or boundaries only ), (not in parallel!!)
 *  ----------------------------------------------------------------------*/


#include "fd.h"
void smooth(int nx,int ny,int nz,float ***mod, float ***mod1){
  
  int i,j,k,ii,jj,kk,l1,l2,l3,iii,jjj,kkk;
  
  int win=20;
  extern int NXG, NYG, NZG,FW, FREE_SURF;
  int smoothpart=0, dum=0, dist=0;
  float a=0.05, dummy1=0.0,dummy2=0.0;

  
	for (k=1;k<=NZG;k++){
		for (i=1;i<=NXG;i++){
			for (j=1;j<=1;j++){
			  
				for(l1=-win;l1<=win;l1++){
				  
					kk=NZG-abs(NZG-abs(k-1+l1)-1);
				
					for(l2=-win;l2<=win;l2++){
					  
						ii=NXG-abs(NXG-abs(i-1+l2)-1);
					
						for(l3=-win;l3<=win;l3++){
							jj=NYG-abs(NYG-abs(j-1+l3)-1);
							
						mod1[j][i][k]+=mod[jj][ii][kk];
						}
					}
				}
			}
			
			for (j=2;j<=NYG;j++){
				 
				mod1[j][i][k]=mod1[j-1][i][k];
			  
				for(l1=-win;l1<=win;l1++){
					kk=NZG-abs(NZG-abs(k-1+l1)-1);
					for(l2=-win;l2<=win;l2++){
						ii=NXG-abs(NXG-abs(i-1+l2)-1);
				  
						jj=NYG-abs(NYG-abs(j-1+win)-1);				
						mod1[j][i][k]+=mod[jj][ii][kk];
				
						jj=NYG-abs(NYG-abs(j-1-1-win)-1);
						mod1[j][i][k]+=-mod[jj][ii][kk];
					}
				}
				mod1[j-1][i][k]=mod1[j-1][i][k]/pow((2*win+1),3);	
				
			}mod1[NYG][i][k]=mod1[NYG][i][k]/pow((2*win+1),3);
		}
	}
	
	/*smooth only boundary*/
	if(smoothpart==1){	
	  
	 /* for (k=FW+1;k<=NZG-FW;k++){
			for (i=FW+1;i<=NXG-FW;i++){
				for (j=dum+1;j<=NYG-FW;j++){
								
					mod1[j][i][k]=mod[j][i][k];
				}
			}
	  }*/
	  
	  dist=FW-2;
		if(FREE_SURF==0)dum=dist;
		for (k=dist+1;k<=NZG-dist;k++){
			kk=k-dist; kkk=NZG-dist+1-k;
			for (i=dist+1;i<=NXG-dist;i++){
				ii=i-dist; iii=NXG-dist+1-i;
				for (j=dum+1;j<=NYG-dist;j++){
					jj=j-dist;jjj=NXG-dist+1-j;
					
					mod1[j][i][k]=mod[j][i][k];
				
					/*mod1[j][i][k]=mod1[dist][i][k]*exp(-a*(jj)*(jj))+mod1[j][i][k]*(1-exp(-a*(jj)*(jj)));*/
					dummy1=mod1[j][dist][k]*exp(-a*(ii)*(ii))+mod1[j][i][k]*(1-exp(-a*(ii)*(ii)));
					dummy2=mod1[j][i][dist]*exp(-a*(kk)*(kk))+dummy1*(1-exp(-a*(kk)*(kk)));
					
					dummy1=mod1[NYG-dist+1][i][k]*exp(-a*(jjj)*(jjj))+dummy2*(1-exp(-a*(jjj)*(jjj)));
					dummy2=mod1[j][NXG-dist+1][k]*exp(-a*(iii)*(iii))+dummy1*(1-exp(-a*(iii)*(iii)));
					dummy1=mod1[j][i][NZG-dist+1]*exp(-a*(kkk)*(kkk))+dummy2*(1-exp(-a*(kkk)*(kkk)));
					mod1[j][i][k]=dummy1;
				  
				}
			}
		}
		
		/*for (k=1;k<=NZG;k++){
			for (i=FW+1;i<=NXG;i++){ii=i-FW;
				for (j=1;j<=NYG;j++){
				  
				  mod1[j][i][k]=mod1[j][FW][k]*exp(-a*(ii)*(ii))+mod1[j][i][k]*(1-exp(-a*(ii)*(ii)));
				}
			}
		}*/
		
	}
	
}


