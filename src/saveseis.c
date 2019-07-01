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
 *   write seismograms to files 
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void saveseis(FILE *fp, float **sectionvx, float **sectionvy,float **sectionvz,
float **sectionp, float **sectioncurl, float **sectiondiv,
int  **recpos, int  **recpos_loc, int ntr, float ** srcpos, int ishot,int ns, int obs, int iteration){ 
		
	extern int SEISMO, SEIS_FORMAT, MYID, RUN_MULTIPLE_SHOTS,NPROC,METHOD;	
	extern char  SEIS_FILE[STRING_SIZE], SEIS_OBS_FILE[STRING_SIZE];
	extern FILE *FP;

	char vxf[STRING_SIZE], vyf[STRING_SIZE], vzf[STRING_SIZE], curlf[STRING_SIZE], divf[STRING_SIZE], pf[STRING_SIZE],file_ext[5];
	char outfile[STRING_SIZE],infile[STRING_SIZE];
	char seisfile[STRING_SIZE];
	float * buf2,* buf1, **srcpos1=NULL;
	int nsrc=1,i,k,m,l,nt;
	FILE *fin, *fout;
	
		srcpos1=fmatrix(1,7,1,1);
		for (nt=1;nt<=7;nt++) srcpos1[nt][1]=srcpos[nt][ishot];
		
		switch (SEIS_FORMAT){
		case 0: sprintf(file_ext,"su"); break;
		case 1: sprintf(file_ext,"su");  break;
		case 2: sprintf(file_ext,"txt"); break;
		case 3: sprintf(file_ext,"bin"); break;
		}
		
		if(obs==0) sprintf(seisfile,"%s",SEIS_FILE);
		if(obs==1) sprintf(seisfile,"%s",SEIS_OBS_FILE);
		
	if(ntr>0){	
		if (RUN_MULTIPLE_SHOTS){
			sprintf(vxf,"%s_vx_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,MYID);
			sprintf(vyf,"%s_vy_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,MYID);
			sprintf(vzf,"%s_vz_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,MYID);
			sprintf(curlf,"%s_rot_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,MYID);
			sprintf(divf,"%s_div_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,MYID);
			sprintf(pf,"%s_p_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,MYID);
		}
		else{
			sprintf(vxf,"%s_vx_it%d.%s.%d",seisfile,iteration,file_ext,MYID);
			sprintf(vyf,"%s_vy_it%d.%s.%d",seisfile,iteration,file_ext,MYID);
			sprintf(vzf,"%s_vz_it%d.%s.%d",seisfile,iteration,file_ext,MYID);
			sprintf(curlf,"%s_rot_it%d.%s.%d",seisfile,iteration,file_ext,MYID);
			sprintf(divf,"%s_div_it%d.%s.%d",seisfile,iteration,file_ext,MYID);
			sprintf(pf,"%s_p_it%d.%s.%d",seisfile,iteration,file_ext,MYID);
			
		}
		

		switch (SEISMO){
		case 1 : /* particle velocities only */
			
			fprintf(fp,"\n PE %d is writing %d seismogramtraces (vx)   to \t %s \n",MYID,ntr,vxf);
			outseis(fp,fopen(vxf,"w"),1,sectionvx,recpos,recpos_loc, ntr,srcpos1,nsrc,ns,SEIS_FORMAT);
			fprintf(fp," PE %d is writing %d seismogramtraces (vy)   to \t %s \n",MYID,ntr,vyf);
			outseis(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc, ntr,srcpos1,nsrc,ns,SEIS_FORMAT);
			fprintf(fp," PE %d is writing %d seismogramtraces (vz)   to \t %s \n",MYID,ntr,vzf);
			outseis(fp,fopen(vzf,"w"),3,sectionvz,recpos,recpos_loc, ntr,srcpos1,nsrc,ns,SEIS_FORMAT);
			
			break;
		case 2 : /* pressure only */
			fprintf(fp," PE %d is writing %d seismogramtraces (p)   to \t %s \n",MYID,ntr,pf);
			outseis(fp,fopen(pf,"w"), 0,sectionp, recpos,recpos_loc, ntr,srcpos1,nsrc,ns,SEIS_FORMAT);
			
			break;
		case 3 : /* curl and div only */
			
			fprintf(fp," PE %d is writing %d seismogramtraces (div)   to \t %s \n",MYID,ntr,divf);
			outseis(fp,fopen(divf,"w"),0,sectiondiv,recpos,recpos_loc,ntr,srcpos1,nsrc,ns,SEIS_FORMAT);
			fprintf(fp," PE %d is writing %d seismogramtraces (curl)  to \t %s \n",MYID,ntr,curlf);
			outseis(fp,fopen(curlf,"w"),0,sectioncurl,recpos,recpos_loc,ntr,srcpos1,nsrc,ns,SEIS_FORMAT);	
			
			break;	
		case 4 : /* everything */
			fprintf(fp,"\n PE %d is writing %d seismogramtraces (vx)   to \t %s \n",MYID,ntr,vxf);
			outseis(fp,fopen(vxf,"w"),1,sectionvx,recpos,recpos_loc, ntr,srcpos1,nsrc,ns,SEIS_FORMAT);
			fprintf(fp," PE %d is writing %d seismogramtraces (vy)   to \t %s \n",MYID,ntr,vyf);
			outseis(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc, ntr,srcpos1,nsrc,ns,SEIS_FORMAT);
			fprintf(fp," PE %d is writing %d seismogramtraces (vz)   to \t %s \n",MYID,ntr,vzf);
			outseis(fp,fopen(vzf,"w"),3,sectionvz,recpos,recpos_loc, ntr,srcpos1,nsrc,ns,SEIS_FORMAT);
			
			fprintf(fp," PE %d is writing %d seismogramtraces (p)    to \t %s \n",MYID,ntr,pf);
			outseis(fp,fopen(pf,"w"), 0,sectionp, recpos,recpos_loc, ntr,srcpos1,nsrc,ns,SEIS_FORMAT);

			fprintf(fp," PE %d is writing %d seismogramtraces (div)  to \t %s \n",MYID,ntr,divf);
			outseis(fp,fopen(divf,"w"),0,sectiondiv,recpos,recpos_loc,ntr,srcpos1,nsrc,ns,SEIS_FORMAT);
			fprintf(fp," PE %d is writing %d seismogramtraces (curl) to \t %s \n",MYID,ntr,curlf);
			outseis(fp,fopen(curlf,"w"),0,sectioncurl,recpos,recpos_loc,ntr,srcpos1,nsrc,ns,SEIS_FORMAT);	
			break;
			
	      }     	
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(MYID==0){
		sprintf(outfile,"%s_vx_it%d.%s.shot%d",seisfile,iteration,file_ext,ishot);
		/*fprintf(FP,"outfile:%s",outfile);*/
		fout=fopen(outfile,"w");
		buf2=vector(1,ns);
		buf1=vector(1,240/sizeof(float));	
			for(i=0;i<=NPROC-1;i++){
				k=0;
				sprintf(infile,"%s_vx_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,i);
				fin=fopen(infile,"r");
				if(fin!=0){
					fseek(fin,0,SEEK_END);
					k=ftell(fin);
					if(k%(240+4*ns)==0)l=k/(240+4*ns);
					else{ fprintf(FP,"attention: Error in seisfile");
						l=0;
					}
					/*fprintf(FP,"l=%d",l);
					fprintf(FP,"infile:%s,MYID=%d\n",infile,i);*/
					fseek(fin,0,SEEK_SET);
					
					for(m=1;m<=l;m++){
					  /*fprintf(FP,"m=%d",m);*/
					  
					  
					fread(buf1,sizeof(float),240/sizeof(float),fin);
					fread(buf2,sizeof(float),ns,fin);
							
					fwrite(buf1,sizeof(float),240/sizeof(float),fout);
					fwrite(buf2,sizeof(float),ns,fout);
					}
					fclose(fin);
					
				}
			}
		fclose(fout);
		free_vector(buf2,1,240/sizeof(float));
		free_vector(buf1,1,ns);
	}

	if(MYID==0){
		sprintf(outfile,"%s_vy_it%d.%s.shot%d",seisfile,iteration,file_ext,ishot);
		/*fprintf(FP,"outfile:%s",outfile);*/
		fout=fopen(outfile,"w");
		buf2=vector(1,ns);
		buf1=vector(1,240/sizeof(float));
		k=0;
		
		for(i=0;i<=NPROC-1;i++){
			sprintf(infile,"%s_vy_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,i);
			fin=fopen(infile,"r");
			if(fin!=0){
				fseek(fin,0,SEEK_END);
				k=ftell(fin);
				if(k%(240+4*ns)==0)l=k/(240+4*ns);
				else{ fprintf(FP,"attention: Error in seisfile");
					l=0;
				}
				/*fprintf(FP,"l=%d",l);
				fprintf(FP,"infile:%s,MYID=%d",infile,MYID);*/
				fseek(fin,0,SEEK_SET);
				
				for(m=1;m<=l;m++){
				  /*fprintf(FP,"m=%d",m);*/
				  
				  
				fread(buf1,sizeof(float),240/sizeof(float),fin);
				fread(buf2,sizeof(float),ns,fin);
						
				fwrite(buf1,sizeof(float),240/sizeof(float),fout);
				fwrite(buf2,sizeof(float),ns,fout);
				}
				fclose(fin);
			}
		}
		fclose(fout);
		free_vector(buf2,1,240/sizeof(float));
		free_vector(buf1,1,ns);
	}
	
	if(MYID==0){
		sprintf(outfile,"%s_vz_it%d.%s.shot%d",seisfile,iteration,file_ext,ishot);
		/*fprintf(FP,"outfile:%s",outfile);*/
		fout=fopen(outfile,"w");
		buf2=vector(1,ns);
		buf1=vector(1,240/sizeof(float));
		k=0;
		for(i=0;i<=NPROC-1;i++){
			sprintf(infile,"%s_vz_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,i);
			fin=fopen(infile,"r");
			if(fin!=0){
				fseek(fin,0,SEEK_END);
				k=ftell(fin);
				if(k%(240+4*ns)==0)l=k/(240+4*ns);
				else{ fprintf(FP,"attention: Error in seisfile");
					l=0;
				}
				/*fprintf(FP,"l=%d",l);
				fprintf(FP,"infile:%s,MYID=%d",infile,MYID);*/
				fseek(fin,0,SEEK_SET);
				for(m=1;m<=l;m++){
				  				  
				  
				fread(buf1,sizeof(float),240/sizeof(float),fin);
				fread(buf2,sizeof(float),ns,fin);
						
				fwrite(buf1,sizeof(float),240/sizeof(float),fout);
				fwrite(buf2,sizeof(float),ns,fout);
				}
				fclose(fin);
				
			}
		}
		fclose(fout);
		free_vector(buf2,1,240/sizeof(float));
		free_vector(buf1,1,ns);
	}
	free_matrix(srcpos1,1,7,1,1);
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(iteration>0&&ntr>0&&METHOD&&!obs){
		sprintf(infile,"%s_vx_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,MYID);
		remove(infile);
		sprintf(infile,"%s_vy_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,MYID);
		remove(infile);
		sprintf(infile,"%s_vz_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,MYID);
		remove(infile);
	}
}
