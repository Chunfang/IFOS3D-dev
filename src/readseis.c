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

/*----------------------------------------------------------------------------
 * reads seismograms from files 
 -----------------------------------------------------------------------------*/


#include "fd.h"

void readseis(int ishot, float **section, float **sectionf, int ntr, int ns, int comp){
  

    extern FILE *FP;
    extern int MYID, RUN_MULTIPLE_SHOTS;
    extern char  SEIS_OBS_FILE[STRING_SIZE];
    char file_ext[5];
    FILE *fpdata;
    char data[STRING_SIZE];
    int i,j;
    
    memset(data, '\0', sizeof(data));    
    //fprintf(FP,"comp=%i",comp);
    
    sprintf(file_ext,"su");

    if (RUN_MULTIPLE_SHOTS){
		if(comp==1){sprintf(data,"%s_vx_it1.%s.shot%d.%d",SEIS_OBS_FILE,file_ext,ishot,MYID);}
		if(comp==2){sprintf(data,"%s_vy_it1.%s.shot%d.%d",SEIS_OBS_FILE,file_ext,ishot,MYID);}
		if(comp==3){sprintf(data,"%s_vz_it1.%s.shot%d.%d",SEIS_OBS_FILE,file_ext,ishot,MYID);}
	}
    else{
	      
		if(comp==1){sprintf(data,"%s_vx.%s.shot%d.%d",SEIS_OBS_FILE,file_ext,ishot,MYID);}
		if(comp==2){sprintf(data,"%s_vy.%s.shot%d.%d",SEIS_OBS_FILE,file_ext,ishot,MYID);}
		if(comp==3){sprintf(data,"%s_vz.%s.shot%d.%d",SEIS_OBS_FILE,file_ext,ishot,MYID);}
	}
    fprintf(FP,"\n reads %s",data);
    
    if(ntr>0){
    fpdata=fopen(data,"r");
     
    for(i=1;i<=ntr;i++){
    
	fseek(fpdata,240,SEEK_CUR);
	
    	for(j=1;j<=ns;j++){
	  section[i][j]=0.0; 
	  fread(&section[i][j],sizeof(float),1,fpdata); }
    }
    
    fclose(fpdata);
    }
    //fprintf(FP," readseis finished\n");

  
}