/*
 *   homogeneous half space 
 *   last update 02.11.02, T. Bohlen
 */

#include "fd.h"

void model(float  ***  rho, float ***  pi, float ***  u, 
float ***  taus, float ***  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, *FL, TAU;
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
	extern char  MFILE[STRING_SIZE];
	extern FILE *FP;
	/* local variables */
	float muv, piv, ws;
	float Vp, Vs, Rho;
	float *pts=NULL, ts=0.0, tp=0.0, sumu, sumpi;
	int i, j, k, l, ii, jj, kk;
	char modfile[STRING_SIZE];
	int checker=10, t1, t2, t3;

	/* parameters for layer 1 */
	const float vp1=5985.0, vs1=3325.0, rho1=2800.0;
      const float vp2=6615.0, vs2=3675.0, rho2=3100.0;

	/*-----------------------------------------------------------------------*/
	sumu=0.0; 
	sumpi=0.0;

	 if(L){
		/* vector for maxwellbodies */
		pts=vector(1,L);
		for (l=1;l<=L;l++) {
			pts[l]=1.0/(2.0*PI*FL[l]);
			eta[l]=DT/pts[l];
		}
		ts=TAU;  
		tp=TAU;
		ws=2.0*PI*FL[1];
		for (l=1;l<=L;l++){
			sumu=sumu+((ws*ws*pts[l]*pts[l]*ts)/(1.0+ws*ws*pts[l]*pts[l]));
			sumpi=sumpi+((ws*ws*pts[l]*pts[l]*tp)/(1.0+ws*ws*pts[l]*pts[l]));
		}
	 }
	
	t1=1;
	for(t1=1;t1<=NZG/checker+1;t1++){
		for(k=(t1-1)*checker+1;k<=t1*checker && k<=NZG ;k++){
			t2=1;
			for(t2=1;t2<=NXG/checker+1;t2++){
				 for(i=(t2-1)*checker+1;i<=t2*checker && i<=NXG ;i++){
					t3=1;
					for(t3=1;t3<=NYG/checker+1;t3++){
						for(j=(t3-1)*checker+1;j<=t3*checker && j<=NYG; j++){
						  
							if((t1+t2+t3)%2!=0){Vp=vp1; Vs=vs1; Rho=rho1;}
							else {Vp=vp2; Vs=vs2; Rho=rho2;}
						  
							muv=Vs*Vs*Rho/(1.0+sumu);
							piv=Vp*Vp*Rho/(1.0+sumpi);

							if ((POS[1]==((i-1)/NX)) && 
							    (POS[2]==((j-1)/NY)) && 
							    (POS[3]==((k-1)/NZ))){
								ii=i-POS[1]*NX;
								jj=j-POS[2]*NY;
								kk=k-POS[3]*NZ;
								if(L){
									taus[jj][ii][kk]=ts;
									taup[jj][ii][kk]=tp;
								}
								u[jj][ii][kk]=muv;
								rho[jj][ii][kk]=Rho;
								pi[jj][ii][kk]=piv;
							}
						}
					}
				 }
			}
		}
	}
	
		  
		
		
	


	
	sprintf(modfile,"%s.mod",MFILE);
	writemod(modfile,rho,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);
         if(L)	free_vector(pts,1,L);
}



