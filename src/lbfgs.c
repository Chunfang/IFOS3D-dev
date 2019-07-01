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

/*-------------------------------------------------------------------------
 * calculatipon of L-BFGS update
 * S. Butzer 2014
 --------------------------------------------------------------------------*/

#include "fd.h"
void lbfgs(float ***grad1, float ***grad2, float ***grad3, float *bfgsscale, float **bfgsmod, float **bfgsgrad, int iteration){

	int m=0,v=0,w=0;
	int i,j,k,l;
	float dummy[2], *dummy1, dummy2=0.0, buf[2]; 
	float *q;
	float h0;
	extern int NX,NY,NZ; /*,HESS*/
	extern FILE *FP;
	extern int BFGSNUM, NUMPAR;

	
	dummy1 = vector(1,BFGSNUM);
	q = vector(1,NUMPAR*NX*NY*NZ);
	dummy[0]=0.0; dummy[1]=0.0;
	buf[0]=0.0; buf[1]=0.0;
		
	m=iteration-BFGSNUM;
	if(m<1) m=1;
	
	fprintf(FP,"Start calculation L-BFGS update");
	
	/*save gradient for bfgsgrad(iteration-1) and set q=gradient, attention grad is -gradient of misfit function*/
	w=(iteration-1)%BFGSNUM;
	if(w==0) w=BFGSNUM;
	l=0;
	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
			for (k=1;k<=NZ;k++){
				l++;
				bfgsgrad[w][l]+=-grad1[j][i][k];
				q[l]=grad1[j][i][k];
			}
		}
	}
	if(NUMPAR==2){
		l=NX*NY*NZ;
		for (j=1;j<=NY;j++){
			for (i=1;i<=NX;i++){
				for (k=1;k<=NZ;k++){
					l++;
					bfgsgrad[w][l]+=-grad2[j][i][k];
					q[l]=grad2[j][i][k];
				}
			}
		}
	}
	
	
	/*calculate bfgsscale and H_0*/
	w=(iteration-1)%BFGSNUM; if(w==0) w=BFGSNUM;
	for (l=1;l<=NUMPAR*NX*NY*NZ;l++){
		dummy[0]+=bfgsgrad[w][l]*bfgsmod[w][l];
		dummy[1]+=bfgsgrad[w][l]*bfgsgrad[w][l];
	}
	
	MPI_Allreduce(&dummy,&buf,2,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
		
	bfgsscale[w]=1/buf[0];
	/*printf("bfgsscale(%d)=%e\n",w,bfgsscale[w]);*/
	h0=buf[0]/buf[1];
	/*printf("h0=%e\n",h0);*/
	/*---------------------------------------*/
	
		
	/*L-BFGS two-loop recursion (numerical optimisation Algorithm 9.1)*/
		  
	l=0;			  
	for(v=iteration-1; v>=m;v--){
	 /* printf("dummy10[%d]=%e",w,dummy1[w]);
	  fprintf(FP,"v1=%d\n",v);*/
		 w=v%BFGSNUM;
		 if(w==0) w=BFGSNUM;
		 fprintf(FP,"w1=%d\n",w);
		 for (l=1;l<=NUMPAR*NX*NY*NZ;l++){
			dummy1[w]+=bfgsscale[w]*bfgsmod[w][l]*q[l];
			/*if(iteration==3)fprintf(FP,"bfgsgrad[%d][%d]=%e,bfgsscale[%d]=%e,bfgsmod[%d][%d]=%e,q[%d]=%e\n", w,l,bfgsgrad[w][l],w,bfgsscale[w],w,l,bfgsmod[w][l],l,q[l]);*/
		 }
		/* printf("dummy11[%d]=%e",w,dummy1[w]);*/
		 buf[0]=0.0;
		 MPI_Allreduce(&dummy1[w],&buf[0],1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
		 dummy1[w]=buf[0];
		 /*printf("dummy1(%d)=%e\n",w,dummy1[w]);*/
		 for (l=1;l<=NUMPAR*NX*NY*NZ;l++){
			q[l]+=-dummy1[w]*bfgsgrad[w][l];
			
		 }
	}
	
	
		 for (l=1;l<=NUMPAR*NX*NY*NZ;l++){
			dummy2=q[l];
			q[l]=dummy2*h0;
		 }
	
	/* if(HESS){l=0;
	   	for (j=1;j<=NY;j++){
			for (i=1;i<=NX;i++){
				for (k=1;k<=NZ;k++){ 
					l++;
					dummy2=q[l];
					q[l]=dummy2*hess[j][i][k];
				}
			}
		}
	 }*/
	 
	for(v=m; v<=iteration-1;v++){
	  /*fprintf(FP,"v2=%d\n",v);*/
		w=v%BFGSNUM;
		if(w==0) w=BFGSNUM;
		/*fprintf(FP,"w2=%d\n",w);*/
		dummy2=0.0;
		
		for (l=1;l<=NUMPAR*NX*NY*NZ;l++){
			dummy2+=bfgsscale[w]*bfgsgrad[w][l]*q[l];
		}
		
		buf[0]=0.0;
		MPI_Allreduce(&dummy2,&buf[0],1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
		dummy2=buf[0];
		
		for (l=1;l<=NUMPAR*NX*NY*NZ;l++){
			q[l]+=bfgsmod[w][l]*(dummy1[w]-dummy2);	
		}
	}
	
	/*save gradients of this iteration and update gradient*/
	w=iteration%BFGSNUM;
		if(w==0) w=BFGSNUM;
	
	l=0;
	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
			for (k=1;k<=NZ;k++){  
				l++;  /*w=(j-1)*Nx*Nz+(i-1)*Nz+k*/  
				bfgsgrad[w][l]=grad1[j][i][k];
				grad1[j][i][k]=q[l];
			}
		}
	}
	if(NUMPAR==2){
		l=NX*NY*NZ;
		for (j=1;j<=NY;j++){
			for (i=1;i<=NX;i++){
				for (k=1;k<=NZ;k++){  
					l++;  /*w=(j-1)*Nx*Nz+(i-1)*Nz+k*/  
					bfgsgrad[w][l]=grad2[j][i][k];
					grad2[j][i][k]=q[l];
				}
			}
		}
	}
free_vector(q,1,NUMPAR*NX*NY*NZ);
free_vector(dummy1,1,m);
}
