/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
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
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
 -----------------------------------------------------------------------------------------*/

/* ----------------------------------------------------------------------
 * This is program IFOS3D.
 * 3D elastic Inversion of Full Observed Seismograms
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "globvar.h"

int main(int argc, char **argv){
	int ns, nt, nseismograms=0, nf1, nf2;
	int lsnap, nsnap=0, lsamp=0, nlsamp=0, buffsize;
	int ntr=0, ntr_loc=0, ntr_glob=0, nsrc=0, nsrc_loc=0, ishot, ishot1, nshots; 

	double 	time1=0.0, time2=0.0, time4=0.0;
	double * time_v_update, * time_s_update, * time_s_exchange,* time_v_exchange, * time_timestep;	
	int * xb, * yb, * zb, l,i,j;
	
	float  *** absorb_coeff=NULL;
	float  ***  sxy=NULL, ***  syz=NULL, ***  sxz=NULL;
	float  ***  sxx=NULL, ***  syy=NULL, ***  szz=NULL;
	float  ***  rxy=NULL, ***  ryz=NULL, ***  rxz=NULL;
	float  ***  rxx=NULL, ***  ryy=NULL, ***  rzz=NULL;
	float  ***  vx=NULL, ***  vy=NULL, ***  vz=NULL;
	float  ***  rho=NULL, ***  pi=NULL, ***  u=NULL;
	float  ***  taus=NULL, ***  taup=NULL, *eta=NULL;
	float  *** uipjp=NULL, *** ujpkp=NULL, *** uipkp=NULL, *** tausipjp=NULL, *** tausjpkp=NULL, *** tausipkp=NULL,*** rjp=NULL, *** rkp=NULL, *** rip=NULL;
	
	float  ** srcpos=NULL, **srcpos_loc=NULL, **srcpos_loc_back=NULL,** srcpos1=NULL, ** signals=NULL;
	int    ** recpos=NULL, ** recpos_loc=NULL, *snum_loc=NULL,*rnum_loc=NULL ;
	float  ** sectionvx=NULL, ** sectionvy=NULL, ** sectionvz=NULL, ** sectionp=NULL,
		** sectioncurl=NULL, ** sectiondiv=NULL;
	
	float *** bufferlef_to_rig, *** bufferrig_to_lef;
	float *** buffertop_to_bot, *** bufferbot_to_top;
	float *** bufferfro_to_bac, *** bufferbac_to_fro;
	
	float *** sbufferlef_to_rig, *** sbufferrig_to_lef;
	float *** sbuffertop_to_bot, *** sbufferbot_to_top;
	float *** sbufferfro_to_bac, *** sbufferbac_to_fro;

	float * K_x=NULL, * alpha_prime_x=NULL, * a_x=NULL, * b_x=NULL, * K_x_half=NULL, * alpha_prime_x_half=NULL, * a_x_half=NULL, * b_x_half=NULL, * K_y=NULL, * alpha_prime_y=NULL, * a_y=NULL, * b_y=NULL, * K_y_half=NULL, * alpha_prime_y_half=NULL, * a_y_half=NULL, * b_y_half=NULL, * K_z=NULL, * alpha_prime_z=NULL, * a_z=NULL, * b_z=NULL, * K_z_half=NULL, * alpha_prime_z_half=NULL, * a_z_half=NULL, * b_z_half=NULL;
	float *** psi_sxx_x=NULL, *** psi_syy_y=NULL, *** psi_szz_z=NULL, *** psi_sxy_y=NULL, *** psi_sxy_x=NULL, *** psi_sxz_x=NULL, *** psi_sxz_z=NULL, *** psi_syz_y=NULL, *** psi_syz_z=NULL, *** psi_vxx=NULL, *** psi_vyy=NULL, *** psi_vzz=NULL, *** psi_vxy=NULL, *** psi_vxz=NULL, *** psi_vyx=NULL, *** psi_vyz=NULL, *** psi_vzx=NULL, *** psi_vzy=NULL;

	float ** sectionread=NULL, ** sectionreadf=NULL, ** sectionvxdiff=NULL, ** sectionvydiff=NULL,** sectionvzdiff=NULL, * misfit;
	float L2=0.0, L2all=0.0, L2f=0.0;
	
	/*inversion variables*/
	/*float  ****  fvx=NULL, ****  fvy=NULL, ****  fvz=NULL;*/
	float *** grad1=NULL, *** grad2=NULL, *** grad3=NULL;
	float *** gradprior1=NULL, *** gradprior2=NULL, *** gradprior3=NULL;
	float *** gradprior4=NULL, *** gradprior5=NULL, *** gradprior6=NULL;
	float *** hess1=NULL, *** hess2=NULL, *** hess3=NULL;
	float **** Ffvx=NULL, **** Ffvy=NULL, **** Ffvz=NULL, **** Ffivx=NULL, **** Ffivy=NULL, **** Ffivz=NULL; 
	float **** Fbvx=NULL, **** Fbvy=NULL, **** Fbvz=NULL, **** Fbivx=NULL, **** Fbivy=NULL, **** Fbivz=NULL;
	int iteration=0, steptest=0,cdf=0, groupnum=0, nf=0, * itpergroup,it_group=0,itmax=0, ntast=1;
	int ntr_hess=0;
	float * step , *finv, *beta;
	float  ***  testrho, ***  testpi, ***  testu;
	float buf;
	/*MPI_Status status;*/
	float dummy=0.0;
	int hloop=0,pshot=0,pshot1=0, pshot_loc=0;
	
	float **bfgsmod1=NULL,**bfgsgrad1=NULL, *bfgsscale1=NULL;     /* *bfgsscale2,*bfgsscale3;**bfgsmod2,**bfgsmod3,**bfgsgrad2, **bfgsgrad3,*/
	
	
	MPI_Request *req_send, *req_rec, *sreq_send, *sreq_rec;
	/* MPI_Status  *send_statuses, *rec_statuses; */
	
	float memdyn, memmodel, memseismograms, membuffer, memtotal,memcpml=0.0,memdynf=0.0, memgrad=0.0, membfgs=0.0;
	float fac1, fac2,fac3;
	char *buff_addr;// ext[10];
	char buffer[STRING_SIZE], bufferstring[10];
	/*char comp[6];*/
	FILE * fpsrc=NULL;

	/* Initialize MPI environment */
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&NP);
	MPI_Comm_rank(MPI_COMM_WORLD,&MYID);

        setvbuf(stdout, NULL, _IONBF, 0);

	/* initialize clock for estimating runtime of program */	
	if (MYID == 0){
		time1=MPI_Wtime();
		clock();
	}	
	
	/* print program name, version, author etc to stdout*/
	if (MYID == 0) info(stdout);
	FDMPIVERS=33; /* 3D isotropic elastic */


	/* PE 0 is reading the parameters from the input file ifos3d.inp */
	if (MYID == 0){
		 if (argc>1) {
		 	strncpy(FILEINP,argv[1],STRING_SIZE);
			fprintf(stderr," \n Input parameter filename read from command line : %s. \n",FILEINP);
		 	if (strlen(FILEINP)>STRING_SIZE-1) {
		 		fprintf(stderr,"\n IFOS cannot handel pathes with more than %d characters.\n",STRING_SIZE-1);
				fprintf(stderr," Error: IFOS could not read input parameter file name. -> Exit. \n\n");
				return -1;
		 	}		 
		 }
		 else {
		 	strcpy(FILEINP,"ifos3d.inp");
			fprintf(stderr," Caution: input parameter filename set to default 'ifos3d.inp'. \n\n");
		 }
		 FP=fopen(FILEINP,"r");
		 
		if (strstr(FILEINP,".json")) {
		//read json formated input file
			read_par_json(stdout, FILEINP);
			fclose(FP);
		} else {
		//read "old" input file *.inp, might not work in future
			err(" Old Input files (.inp) are no longer supported. \n Please use .json input files instead. \n\n");
		}
		 
	} 
	       
	         
	/* PE 0 will broadcast the parameters to all others PEs */
	exchange_par(); 
			
	/* Print info on log-files to stdout */
	if (MYID == 0) note(stdout);
	
	/* open log-file (each PE is using different file) */
	/*	fp=stdout; */
	/*sprintf(ext,".%i",MYID);  
	strcat(LOG_FILE,ext);*/	

	/* nodes MYIDo writes logging info to LOG_FILE or stdout */	
	if (MYID==0) FP=stdout; /* logging information will be written to standard output */
		
	/* all other nodes write logging info to LOG_FILE */		
	if (MYID>0) {
	/*if ((FP=fopen(LOG_FILE,"w"))==NULL) err(" Opening log-file failed.");
	fprintf(FP," This is the log-file generated by PE %d \n\n",MYID);*/
	FP=fopen("/dev/null","w");
	FI=fopen("/dev/null","w");	}
	
	fprintf(FP," This is the log-file generated by PE %d \n\n",MYID);
	
	if (MYID==0){
		FI=fopen("in_and_out/ifos3D_invers.out","w");
		setvbuf(FI, NULL, _IONBF, 0);
		fprintf(FI,"------------------ Inversion Parameters ----------------- \n");
	}

	/* domain decomposition */
	initproc();

	/* set some time counters */
	NT=(int)ceil(TIME/DT); /* number of timesteps - replaces: NT=iround(TIME/DT); */
	TIME=(NT-1)*DT; /* TIME set to true time of the last time step */
	if (NDTSHIFT>NT) ns=0;
	else ns=(int)ceil((float)(NT-NDTSHIFT)/(float)NDT); /* number of samples per trace - replaces buggy formula: ns=iround(NT-NDTSHIFT/NDT); */
	lsnap=iround(TSNAP1/DT); /* first snapshot at this timestep */
	
	/* output of parameters to stdout: */
	if (MYID==0) writepar(FP,ns);

	/* NXG, NYG NZG denote size of the entire (global) grid */
	NXG=NX;
	NYG=NY;
	NZG=NZ;

	/* In the following, NX, MY, NZ denote size of the local grid ! */
	NX = IENDX;
	NY = IENDY;
	NZ = IENDZ;

	/* compute receiver locations within each subgrid and
	   store local receiver coordinates in recpos_loc */	
	if (SEISMO){
		fprintf(FP,"\n ------------------ READING RECEIVER PARAMETERS ----------------- \n");
		recpos=receiver(FP,&ntr);
		rnum_loc=ivector(1,ntr);
		recpos_loc = splitrec(recpos,&ntr_loc, ntr, rnum_loc);
		ntr_glob=ntr;
		ntr=ntr_loc;
	}
	if (METHOD) srcpos_loc_back=fmatrix(1,7,1,ntr_loc);
// 		

	/* number of seismogram sections which have to be stored in core memory*/
	switch (SEISMO){
	case 1 : /* particle velocities only */
		nseismograms=3;	
	break;	
	case 2 : /* pressure only */
		nseismograms=1;	
	break;	
	case 3 : /* curl and div only */
		nseismograms=2;		
	break;	
	case 4 : /* everything */
		nseismograms=6;		
	break;
}	
if(METHOD) nseismograms+=4;

	/*allocate memory for dynamic, static and buffer arrays */
	fac1=(NZ+FDORDER)*(NY+FDORDER)*(NX+FDORDER);
	fac2=sizeof(float)*pow(2.0,-20.0);
	fac3=NZ*NY*NX;

	 if (L>0){ /*viscoelastic case*/
            memdyn=15.0*fac1*fac2;
            memmodel=16.0*fac1*fac2;
        }
        else { /* elastic case*/
            memdyn=9.0*fac1*fac2;
            memmodel=10.0*fac1*fac2;
        }
	memseismograms=nseismograms*ntr_glob*ns*fac2;
	membuffer=4.0*6.0*((NX*NZ)+(NY*NZ)+(NX*NY))*fac2;
	buffsize=(FDORDER)*4.0*6.0*(max((NX*NZ),max((NY*NZ),(NX*NY))))*sizeof(MPI_FLOAT);
	if (ABS_TYPE==1) memcpml=2.0*FW*6.0*(NY*NZ+NX*NZ+NY*NX)*fac2+24.0*2.0*FW*fac2;
	if(METHOD){
		memgrad=12*fac2*fac3;
		memdynf=12*NFMAX*(ntr_hess+1)*fac1*fac2;
		if(HESS) memgrad+=3*fac2*fac3;
		if(LBFGS) membfgs=NUMPAR*BFGSNUM*3*fac3;
	}
	memtotal=memdyn+memmodel+memseismograms+membuffer+(buffsize*pow(2.0,-20.0))+memgrad+memdynf+membfgs+memcpml;
	
	if (MYID==0){
		fprintf(FP,"\n ------------------ MEMORY ALLOCATION --------------------------- \n");
		fprintf(FP,"\n **Message from main (printed by PE %d):\n",MYID);
		fprintf(FP," Size of local grids: NX=%d \t NY=%d \t NZ=%d \n",NX,NY,NZ);
		fprintf(FP," Each process is now trying to allocate memory for:\n");
		fprintf(FP," Dynamic modeling variables: \t\t %6.2f MB\n", memdyn);
		fprintf(FP," Static modeling variables: \t\t %6.2f MB\n", memmodel);
   		fprintf(FP," Seismograms: \t\t\t\t %6.2f MB\n", memseismograms);
		fprintf(FP," Dynamic inversion variables: \t\t %6.2f MB\n", memdynf);
		fprintf(FP," Static inversion variables: \t\t %6.2f MB\n", memgrad+membfgs);
		fprintf(FP," Buffer arrays for grid exchange: \t %6.2f MB\n", membuffer);
		fprintf(FP," Network Buffer for MPI_Bsend: \t\t %6.2f MB\n", buffsize*pow(2.0,-20.0));
		fprintf(FP," ------------------------------------------------ \n");
		fprintf(FP," Total memory required: \t\t %6.2f MB.\n\n", memtotal);
	}

MPI_Barrier(MPI_COMM_WORLD); 
	/* allocate buffer for buffering messages */
	buff_addr=malloc(buffsize);
	if (!buff_addr) err("allocation failure for buffer for MPI_Bsend !");
	MPI_Buffer_attach(buff_addr,buffsize);
	
	
	/* allocation for request and status arrays */
	req_send=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
	req_rec=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
	sreq_send=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
	sreq_rec=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
	/* send_statuses=(MPI_Status *)malloc(REQUEST_COUNT*sizeof(MPI_Status));
	rec_statuses=(MPI_Status *)malloc(REQUEST_COUNT*sizeof(MPI_Status)); */
	
        /* allocation for timing arrays used for performance analysis */
	time_v_update=dvector(1,NT);
	time_s_update=dvector(1,NT);
	time_s_exchange=dvector(1,NT);
	time_v_exchange=dvector(1,NT);
	time_timestep=dvector(1,NT);
	
	l=1;
	if(ABS_TYPE==1 && FDORDER==2){l=2;}
	
	if(HESS&&!READ_HESS) ntr_hess=ntr_glob/REC_HESS;
	else     ntr_hess=0;
	/* memory allocation for dynamic (wavefield) arrays */
	if(POS[2]==0){
		vx  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		vy  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		vz  =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		/*fvx  =  f4tensor(1,NT,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		fvy  =  f4tensor(1,NT,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		fvz  =  f4tensor(1,NT,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);*/
		if(METHOD){
			Ffvx  =  f4tensor(1,NFMAX,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Ffvy  =  f4tensor(1,NFMAX,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Ffvz  =  f4tensor(1,NFMAX,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Fbvx  =  f4tensor(1,NFMAX*(ntr_hess+1),0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Fbvy  =  f4tensor(1,NFMAX*(ntr_hess+1),0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Fbvz  =  f4tensor(1,NFMAX*(ntr_hess+1),0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Ffivx  =  f4tensor(1,NFMAX,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Ffivy  =  f4tensor(1,NFMAX,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Ffivz  =  f4tensor(1,NFMAX,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Fbivx  =  f4tensor(1,NFMAX*(ntr_hess+1),0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Fbivy  =  f4tensor(1,NFMAX*(ntr_hess+1),0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Fbivz  =  f4tensor(1,NFMAX*(ntr_hess+1),0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		}
		
		sxy =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		syz =  f3tensor(0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	}

	if(POS[2]>0){
		vx  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		vy  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		vz  =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		/*fvx  =  f4tensor(1,NT,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		fvy  =  f4tensor(1,NT,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		fvz  =  f4tensor(1,NT,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);*/
		if(METHOD){
			Ffvx  =  f4tensor(1,NFMAX,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Ffvy  =  f4tensor(1,NFMAX,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Ffvz  =  f4tensor(1,NFMAX,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Fbvx  =  f4tensor(1,NFMAX*(ntr_hess+1),1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Fbvy  =  f4tensor(1,NFMAX*(ntr_hess+1),1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Fbvz  =  f4tensor(1,NFMAX*(ntr_hess+1),1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Ffivx  =  f4tensor(1,NFMAX,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Ffivy  =  f4tensor(1,NFMAX,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Ffivz  =  f4tensor(1,NFMAX,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Fbivx  =  f4tensor(1,NFMAX*(ntr_hess+1),1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Fbivy  =  f4tensor(1,NFMAX*(ntr_hess+1),1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
			Fbivz  =  f4tensor(1,NFMAX*(ntr_hess+1),1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		}
		sxy =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		syz =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	}

	sxz =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	sxx =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	syy =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	szz =  f3tensor(1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	
	xb=ivector(0,1);
	yb=ivector(0,1);
	zb=ivector(0,1);

	/*memory allocation for inversion parameter*/
	if(METHOD){
		grad1 = f3tensor(1,NY,1,NX,1,NZ);
		grad2 = f3tensor(1,NY,1,NX,1,NZ);
		grad3 = f3tensor(1,NY,1,NX,1,NZ);
		gradprior1 = f3tensor(1,NY,1,NX,1,NZ);
		gradprior2 = f3tensor(1,NY,1,NX,1,NZ);
		gradprior3 = f3tensor(1,NY,1,NX,1,NZ);
		gradprior4 = f3tensor(1,NY,1,NX,1,NZ);
		gradprior5 = f3tensor(1,NY,1,NX,1,NZ);
		gradprior6 = f3tensor(1,NY,1,NX,1,NZ);
		if(HESS){
			hess1 = f3tensor(1,NY,1,NX,1,NZ);
			hess2 = f3tensor(1,NY,1,NX,1,NZ);
			hess3 = f3tensor(1,NY,1,NX,1,NZ);
		}
	}
	
	if(LBFGS){
		bfgsmod1=fmatrix(1,BFGSNUM,1,NUMPAR*NX*NY*NZ);
		/*bfgsmod2=fmatrix(1,bfgsnum,1,NX*NY*NZ);
		bfgsmod3=fmatrix(1,bfgsnum,1,NX*NY*NZ);*/
		bfgsgrad1=fmatrix(1,BFGSNUM,1,NUMPAR*NX*NY*NZ);
		/*bfgsgrad2=fmatrix(1,bfgsnum,1,NX*NY*NZ);
		bfgsgrad3=fmatrix(1,bfgsnum,1,NX*NY*NZ);*/
		bfgsscale1=vector(1,BFGSNUM);
		/*bfgsscale2=vector(1,bfgsnum);
		bfgsscale3=vector(1,bfgsnum);*/
	}
	
	misfit = vector(0,3);
	step = vector(0,4);
	beta=vector(0,2);
	finv=vector(0,NFMAX-1);
	itpergroup=ivector(0,1);
		MPI_Barrier(MPI_COMM_WORLD);
	/* memory allocation for CPML variables*/
        if(ABS_TYPE==1){

        K_x = vector(1,2*FW);
        alpha_prime_x = vector(1,2*FW);
        a_x = vector(1,2*FW);
        b_x = vector(1,2*FW);
        K_x_half = vector(1,2*FW);
        alpha_prime_x_half = vector(1,2*FW);
        a_x_half = vector(1,2*FW);
        b_x_half = vector(1,2*FW);

        K_y = vector(1,2*FW);
        alpha_prime_y = vector(1,2*FW);
        a_y = vector(1,2*FW);
        b_y = vector(1,2*FW);
        K_y_half = vector(1,2*FW);
        alpha_prime_y_half = vector(1,2*FW);
        a_y_half = vector(1,2*FW);
        b_y_half = vector(1,2*FW);

        K_z = vector(1,2*FW);
        alpha_prime_z = vector(1,2*FW);
        a_z = vector(1,2*FW);
        b_z = vector(1,2*FW);
        K_z_half = vector(1,2*FW);
        alpha_prime_z_half = vector(1,2*FW);
        a_z_half = vector(1,2*FW);
        b_z_half = vector(1,2*FW);


        psi_sxx_x =  f3tensor(1,NY,1,2*FW,1,NZ); 
        psi_sxy_x =  f3tensor(1,NY,1,2*FW,1,NZ);
        psi_sxz_x =  f3tensor(1,NY,1,2*FW,1,NZ);
        psi_syy_y =  f3tensor(1,2*FW,1,NX,1,NZ);
        psi_sxy_y =  f3tensor(1,2*FW,1,NX,1,NZ);
        psi_syz_y =  f3tensor(1,2*FW,1,NX,1,NZ);
        psi_szz_z =  f3tensor(1,NY,1,NX,1,2*FW);
        psi_sxz_z =  f3tensor(1,NY,1,NX,1,2*FW);
        psi_syz_z =  f3tensor(1,NY,1,NX,1,2*FW);


        psi_vxx   =  f3tensor(1,NY,1,2*FW,1,NZ);
        psi_vyy   =  f3tensor(1,2*FW,1,NX,1,NZ);
        psi_vzz   =  f3tensor(1,NY,1,NX,1,2*FW);
        psi_vxy   =  f3tensor(1,2*FW,1,NX,1,NZ);
        psi_vxz   =  f3tensor(1,NY,1,NX,1,2*FW);
        psi_vyx   =  f3tensor(1,NY,1,2*FW,1,NZ);
        psi_vyz   =  f3tensor(1,NY,1,NX,1,2*FW);
        psi_vzx   =  f3tensor(1,NY,1,2*FW,1,NZ);
        psi_vzy   =  f3tensor(1,2*FW,1,NX,1,NZ);

        }

	if(L){
		rxy =  f3tensor(1,NY,1,NX,1,NZ);
		ryz =  f3tensor(1,NY,1,NX,1,NZ);
		rxz =  f3tensor(1,NY,1,NX,1,NZ);
		rxx =  f3tensor(1,NY,1,NX,1,NZ);
		ryy =  f3tensor(1,NY,1,NX,1,NZ);
		rzz =  f3tensor(1,NY,1,NX,1,NZ);
		taus=  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
		taup=  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
		eta =  vector(1,L);
		tausipjp=f3tensor(1,NY,1,NX,1,NZ);
		tausjpkp=f3tensor(1,NY,1,NX,1,NZ);
		tausipkp=f3tensor(1,NY,1,NX,1,NZ);
	}
	
	/* memory allocation for static (model) arrays */
	rho =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
	pi  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
	u   =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
	absorb_coeff=  f3tensor(1,NY,1,NX,1,NZ);
	testrho =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
	testpi  =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
	testu   =  f3tensor(0,NY+1,0,NX+1,0,NZ+1);
	
	/* averaged material parameters */
	uipjp=f3tensor(1,NY,1,NX,1,NZ);
	ujpkp=f3tensor(1,NY,1,NX,1,NZ);
	uipkp=f3tensor(1,NY,1,NX,1,NZ);
	
	rjp=f3tensor(1,NY,1,NX,1,NZ);
	rkp=f3tensor(1,NY,1,NX,1,NZ);
	rip=f3tensor(1,NY,1,NX,1,NZ);
	

	
	/* memory allocation for buffer arrays in which the wavefield
	   information which is exchanged between neighbouring PEs is stored */
	
	/* number of wavefield parameters that need to be exchanged - see exchange_v.c */
	nf1=(3*FDORDER/2)-1;
	nf2=nf1-1;
	
	bufferlef_to_rig = f3tensor(1,NY,1,NZ,1,nf1);
	bufferrig_to_lef = f3tensor(1,NY,1,NZ,1,nf2);
	buffertop_to_bot = f3tensor(1,NX,1,NZ,1,nf1);
	bufferbot_to_top = f3tensor(1,NX,1,NZ,1,nf2);
	bufferfro_to_bac = f3tensor(1,NY,1,NX,1,nf1);
	bufferbac_to_fro = f3tensor(1,NY,1,NX,1,nf2);
	
	sbufferlef_to_rig = f3tensor(1,NY,1,NZ,1,nf2);
	sbufferrig_to_lef = f3tensor(1,NY,1,NZ,1,nf1);
	sbuffertop_to_bot = f3tensor(1,NX,1,NZ,1,nf2);
	sbufferbot_to_top = f3tensor(1,NX,1,NZ,1,nf1);
	sbufferfro_to_bac = f3tensor(1,NY,1,NX,1,nf2);
	sbufferbac_to_fro = f3tensor(1,NY,1,NX,1,nf1);

	if (ntr>0){
	switch (SEISMO){
	case 1 : /* particle velocities only */
		sectionvx=fmatrix(1,ntr,1,ns);
		sectionvy=fmatrix(1,ntr,1,ns);	
		sectionvz=fmatrix(1,ntr,1,ns);	
		break;	
	case 2 : /* pressure only */
		sectionp=fmatrix(1,ntr,1,ns);
		break;	
	case 3 : /* curl and div only */
		sectioncurl=fmatrix(1,ntr,1,ns);
		sectiondiv=fmatrix(1,ntr,1,ns);	
		break;	
	case 4 : /* everything */
		sectionvx=fmatrix(1,ntr,1,ns);
		sectionvy=fmatrix(1,ntr,1,ns);	
		sectionvz=fmatrix(1,ntr,1,ns);	
		sectioncurl=fmatrix(1,ntr,1,ns);
		sectiondiv=fmatrix(1,ntr,1,ns);		
		sectionp=fmatrix(1,ntr,1,ns);
		break;
	}
	}

	/* memory for inversion */
	sectionread=fmatrix(1,ntr_glob,1,ns);
	sectionreadf=fmatrix(1,ntr_glob,1,1);
	sectionvxdiff=fmatrix(1,ntr_glob,1,ns);
	sectionvydiff=fmatrix(1,ntr_glob,1,ns);
	sectionvzdiff=fmatrix(1,ntr_glob,1,ns);
	
        /* memory for source position definition */
	srcpos1=fmatrix(1,6,1,1);
	
	

	if (MYID==0) 
		fprintf(FP," ... memory allocation for PE %d was successfull.\n\n", MYID);

	
	/* create model grids */
	fprintf(FP,"\n ------------------ MODEL CREATION AND OUTPUT-------------------- \n");
	if(READMOD) readmod(rho,pi,u,taus,taup,eta);
	else model(rho,pi,u,taus,taup,eta);	
	/*model_gauss(rho,pi,u,taus,taup,eta);*/
	/*model2_5(rho,pi,u,taus,taup);*/
	outmod(NX,NY,NZ,rho,pi,u,0);
	/*smooth(NX,NY,NZ,pi,testpi);
	smooth(NX,NY,NZ,u,testu);
	smooth(NX,NY,NZ,rho,testrho);*/
		
	if(ABS_TYPE==1){    
	 CPML_ini_elastic(xb,yb,zb);
	 }

	if(ABS_TYPE==2){    
	xb[0]=1; xb[1]=NX;
	yb[0]=1; yb[1]=NY;
	zb[0]=1; zb[1]=NZ;
        }    
        
	
	
	
	/* Reading source positions from SOURCE_FILE */ 	
	fprintf(FP,"\n ------------------ READING SOURCE PARAMETERS ------------------- \n");
	
	    if (MYID==0) switch (SRCREC) {
		case 0: 
			if (MYID==0) err("SRCREC parameter is invalid (SRCREC!=1)! No source parameters specified!");
			break;
		case 1:
			fprintf(FP,"\n Reading source parameters from file: %s (IFOS source format)\n",SOURCE_FILE);
			if ((fpsrc=fopen(SOURCE_FILE,"r"))==NULL) err(" Source file could not be opened !");
			while(fgets(buffer, STRING_SIZE, fpsrc))
			{
			        sscanf(buffer,"%s",bufferstring);
				/* checks if the line contains a '%'character which indicates a comment line,
					and if the reading of a string was successful, which is not the case for an empty line*/
				if ((strchr(buffer,'%')==0) && (sscanf(buffer,"%s",bufferstring)==1)) ++(nsrc);			
   		   	}
			rewind(fpsrc);
			if ((nsrc)==0) fprintf(FP,"\n WARNING: Could not determine number of sources parameter sets in input file. Assuming %d.\n",(nsrc=0));
			else fprintf(FP," Number of source positions specified in %s : %d \n",SOURCE_FILE,nsrc);
			break;
		case 2: if ((PLANE_WAVE_DEPTH>0)) {
				/*determining the number of sources in the specified plane normal/tilted to the surface/upper model boundary*/
				nsrc=(NXG-2*FW+1)*(NYG-2*FW+1);
				/*fprintf(FP,"\n nsrc= %i with NGX=%i, NYG=%i and FW=%i. \n",nsrc,NXG,NYG,FW);*/
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Bcast(&nsrc,1,MPI_INT,0,MPI_COMM_WORLD);
				srcpos= fmatrix(1,7,1,nsrc);
				pwsources(&nsrc,srcpos);	
			} else {
				err("SRCREC parameter specifies PLANE_WAVE excitation, but PLANE_WAVE_DEPTH<=0!");
			}	
			break;
	
		default: err("SRCREC parameter is invalid (SRCREC!=1 or SRCREC!=2)! No source parameters specified!");
           }
	   
	   MPI_Bcast(&nsrc,1,MPI_INT,0,MPI_COMM_WORLD);
	   MPI_Barrier(MPI_COMM_WORLD);
	   srcpos= fmatrix(1,7,1,nsrc); 
	   sources(fpsrc,&nsrc,srcpos);
	   /*originally, SOURCE_TYPE (stype) is defined in the source file, if not, SOURCE_TYPE is taken from the input file */
	   /*if (stype==NULL) printf("PE%d: Source type(s) undefined?! \n",MYID);*/
	
	
	snum_loc = ivector(1,nsrc);
	srcpos_loc = splitsrc(srcpos,&nsrc_loc, nsrc,snum_loc);
	
	if (RUN_MULTIPLE_SHOTS){
		nshots=nsrc; 
		if (nsrc_loc>0) signals=fmatrix(1,1,1,NT);}
	else   {nshots=1;
		if (nsrc_loc>0) signals=fmatrix(1,nsrc_loc,1,NT);}
	
	
	/* check if the FD run will be stable and free of numerical dispersion */
	checkfd(FP,rho,pi,u,taus,taup,eta,srcpos,nsrc,recpos,ntr_glob);
     
	
	/* calculate damping coefficients for CPML boundary*/
          if(ABS_TYPE==1){   
CPML_coeff(K_x,alpha_prime_x,a_x,b_x,K_x_half,alpha_prime_x_half,a_x_half,b_x_half,K_y,alpha_prime_y,a_y,b_y,K_y_half,alpha_prime_y_half,a_y_half,b_y_half,K_z,alpha_prime_z,a_z,b_z,K_z_half,alpha_prime_z_half,a_z_half,b_z_half);
          }
          
	/* calculate 3-D array for exponential damping of reflections
	   at the edges of the numerical mesh */
	if(ABS_TYPE==2){   
	  absorb(absorb_coeff);
        }


	/* comunication initialisation for persistent communication */
	/*comm_ini(bufferlef_to_rig, bufferrig_to_lef,
	    buffertop_to_bot, bufferbot_to_top, bufferfro_to_bac,
	    bufferbac_to_fro, req_send, req_rec);
	    
	comm_ini_s(sbufferlef_to_rig, sbufferrig_to_lef,
	    sbuffertop_to_bot, sbufferbot_to_top, sbufferfro_to_bac,
	    sbufferbac_to_fro, sreq_send, sreq_rec);*/
	 
	 /* initialisation of PML and ABS domain */
	if(ABS_TYPE==1){    
	 CPML_ini_elastic(xb,yb,zb);
	 }

	if(ABS_TYPE==2){    
	xb[0]=1; xb[1]=NX;
	yb[0]=1; yb[1]=NY;
	zb[0]=1; zb[1]=NZ;
        }    
        
	
	if(HESS&&READ_HESS==0) hloop=ntr_hess;
	if(HESS&&READ_HESS==1) readhess(NX,NY,NZ,hess1,hess2,hess3,finv[0], iteration+1-it_group);
	MPI_Barrier(MPI_COMM_WORLD);
	if(METHOD==1 && EXTOBS==1){
		for (ishot=1;ishot<=nshots;ishot++){
		if(ntr_loc>0){
			readseis_split(ishot, sectionvx, ntr_glob, rnum_loc, ns, 1);
			readseis_split(ishot, sectionvy, ntr_glob, rnum_loc, ns, 2);
			readseis_split(ishot, sectionvz, ntr_glob, rnum_loc, ns, 3);}
			saveseis(FP,sectionvx,sectionvy,sectionvz,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos,ishot,ns,1,1);	
		}
	}
	for (iteration=ITMIN;iteration<=ITMAX;iteration++){
	  
		if (MYID==0){
			time2=MPI_Wtime();
			fprintf(FI,"\n *********** STARTING ITERATION NUMBER %d ***************\n",iteration);
			fprintf(FP,"\n *********** STARTING ITERATION NUMBER %d ***************\n",iteration);
		}
		
		if(METHOD){
			  zero_grad(NX,NY,NZ,grad1,grad2,grad3);  
			  zero_invers(NX,NY,NZ,Ffvx,Ffvy,Ffvz,Ffivx,Ffivy,Ffivz,Fbvx,Fbvy,Fbvz,Fbivx,Fbivy,Fbivz,NFMAX,ntr_hess);}
		
		matcopy(rho,pi,u,taus,taup);
		if (FREE_SURF)constant_boundary(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],pi, u, rho);
		av_mat(rho,pi,u,taus,taup,uipjp,ujpkp,uipkp,tausipjp,tausjpkp,tausipkp,rjp,rkp,rip);
		L2all=0.0;
		L2f=0.0;
		misfit[0]=0.0;
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		if(METHOD==1){
			if(it_group>=itmax||iteration==1){
				readinv(finv,&nf,&groupnum,itpergroup,NFMAX);

				cdf=1;
				itmax=itpergroup[1];
				it_group=1;
			}
			else {
				cdf=0;
				it_group+=1;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
				
		if (MYID==0){fprintf(FI,"cdf=%d, f=%4.2f, itmax=%d \n",cdf,finv[0],itmax);}
		
		if(METHOD){
			for(i=1;i<=4;i++){
			for(j=1;j<=ntr_loc;j++){srcpos_loc_back[i][j]=recpos_loc[i][j];	}}
		}
		
		for (ishot=1;ishot<=nshots;ishot++){
			L2=0.0; 
		  
			if (MYID==0) fprintf(FI,"SHOT %d",ishot);
			fprintf(FP,"\n **********   Starting simulation for shot %d of %d  ********** \n",ishot,nshots);
			
			if (RUN_MULTIPLE_SHOTS){nsrc_loc=snum_loc[ishot];
				if(nsrc_loc>0){
					for (nt=4;nt<=7;nt++) srcpos_loc[nt][1]=srcpos[nt][ishot];
					srcpos_loc[1][1]=(float)(((iround(srcpos[1][ishot]/DX)-1)%IENDX)+1);
					srcpos_loc[2][1]=(float)(((iround(srcpos[2][ishot]/DY)-1)%IENDY)+1);
					srcpos_loc[3][1]=(float)(((iround(srcpos[3][ishot]/DZ)-1)%IENDZ)+1);
				}
			}
			/*printf("source=%e,%e,%e,%e,%e,%e,%e\n",srcpos_loc[1][1],srcpos_loc[2][1],srcpos_loc[3][1],srcpos_loc[4][1],srcpos_loc[5][1],srcpos_loc[6][1],srcpos_loc[7][1]);*/
			if(nsrc_loc>0){	wavelet(srcpos_loc,nsrc_loc,SOURCE_SHAPE,signals);	
			/*if (SEISMO){
			  	saveseis(FP,signals,sectionvy,sectionvz,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos,ishot,ns,0,11);
			}*/
			if(FILT==1) filt_seis(signals,nsrc_loc,NT,finv[nf-1]);
			/*if (SEISMO){
			  	saveseis(FP,signals,sectionvy,sectionvz,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos,ishot,ns,0,12);
			}*/
			}
						
			/* initialize wavefield with zero */		
			zero_wavefield(NX,NY,NZ,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rxx,ryy,rzz,rxy,ryz,rxz,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_sxy_y,psi_syy_y,psi_syz_y,psi_sxz_z,psi_syz_z,psi_szz_z,psi_vxx,psi_vyx,psi_vzx,psi_vxy,psi_vyy,psi_vzy,psi_vxz,psi_vyz,psi_vzz);			
			if(METHOD) zero_invers(NX,NY,NZ,Ffvx,Ffvy,Ffvz,Ffivx,Ffivy,Ffivz,Fbvx,Fbvy,Fbvz,Fbivx,Fbivy,Fbivz,NFMAX,0);
			
			/*---------------------------- forward propagation-----------------------------------------------*/
			pshot=0;
			lsamp=NDTSHIFT+1;
			nlsamp=1;
			if(METHOD==1){
				dummy=(1/(finv[nf-1]*TAST*DT));
				ntast=(dummy);
				if(!ntast) ntast=1;
				//fprintf(FP, "ntast=%i, TAST=%i \n", ntast, TAST);
			}
			
			if(MYID==0) fprintf(FP,"\n ****************************************\n ");
			
			for (nt=1;nt<=NT;nt++){
				if(MYID==0) if(!(nt%(NT/40))) fprintf(FP,"*");
				
				time_v_update[nt]=0.0;
				time_s_update[nt]=0.0;

				/* Check if simulation is still stable */
				if (isnan(vy[NY/2][NX/2][NZ/2])) err(" Simulation is unstable !"); /* maybe just breaking the loop would be better */
								
				/* update of particle velocities */	
				time_v_update[nt]+=update_v(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],nt,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rho,rjp, rkp,rip,srcpos_loc,signals,signals,signals,nsrc_loc,absorb_coeff,0);

				if(ABS_TYPE==1){
				update_v_CPML(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],nt,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rho,rjp,rkp,rip,srcpos_loc,signals,signals,signals,nsrc_loc,absorb_coeff,0,K_x,a_x,b_x,K_x_half,a_x_half,b_x_half, K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,K_z,a_z,b_z,K_z_half,a_z_half,b_z_half,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_sxy_y,psi_syy_y,psi_syz_y,psi_sxz_z,psi_syz_z,psi_szz_z);}
		
				/* exchange values of particle velocities at grid boundaries between PEs */
				time_v_exchange[nt]+=exchange_v(vx,vy,vz, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, 
								bufferfro_to_bac, bufferbac_to_fro, req_send, req_rec);

				/*time_s_update[nt]+=update_s(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rxx,ryy,rzz,rxy,ryz,rxz,pi,u,uipjp,ujpkp,uipkp,taus,tausipjp,tausjpkp,tausipkp,taup,eta);*/
				time_s_update[nt]+=update_s_elastic(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,pi,u,uipjp,ujpkp,uipkp);

				if(ABS_TYPE==1){
				/*update_s_CPML(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rxx,ryy,rzz,rxy,ryz,rxz,pi,u,uipjp,ujpkp,uipkp,taus,tausipjp,tausjpkp,tausipkp,taup,eta,K_x,a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,K_z,a_z,b_z,K_z_half,a_z_half,b_z_half,psi_vxx,psi_vyx,psi_vzx,psi_vxy,psi_vyy,psi_vzy,psi_vxz,psi_vyz,psi_vzz);*/
				update_s_CPML_elastic(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,pi,u,uipjp,ujpkp,uipkp,K_x,a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,K_z,a_z,b_z,K_z_half,a_z_half,b_z_half,psi_vxx,psi_vyx,psi_vzx,psi_vxy,psi_vyy,psi_vzy,psi_vxz,psi_vyz,psi_vzz);
				}
				
				/* explosive source */
				psource(nt,sxx,syy,szz,srcpos_loc,signals,nsrc_loc);
				
				/* exchange values of stress at boundaries between PEs */
				time_s_exchange[nt]+=exchange_s(sxx,syy,szz,sxy,syz,sxz,sbufferlef_to_rig,sbufferrig_to_lef,sbuffertop_to_bot,sbufferbot_to_top,sbufferfro_to_bac,sbufferbac_to_fro, sreq_send, sreq_rec);
						
				/* extracting monochromatic wavefields */
				if(nt%ntast==0&&METHOD){
					discfourier(1,NX,1,NY,1,NZ,nt,vx,vy,vz,Ffvx,Ffvy,Ffvz, Ffivx, Ffivy,Ffivz,finv,nf,ntast,pshot,0);
					l++;
				}

				/* stress free surface ? */
				if ((FREE_SURF) && (POS[2]==0))
					surface_elastic(1,u,pi,sxx,syy,szz,sxy,syz,vx,vy,vz);
					/*surface(1,u,pi,taus,taup,eta,sxx,syy,szz,sxy,syz,rxx,ryy,rzz,vx,vy,vz);*/		

				/* store amplitudes at receivers in sectionvx-sectionvz */
				if ((SEISMO) && (ntr>0) && (nt==lsamp)){
					seismo(nlsamp,ntr,recpos_loc,sectionvx,sectionvy,sectionvz,
						sectiondiv,sectioncurl,sectionp,vx,vy,vz,sxx,syy,szz,pi,u);
					nlsamp++;
					lsamp+=NDT;
				}     		
				
				/* save snapshot in file */
				if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){
				    MPI_Barrier(MPI_COMM_WORLD);
				    snap(FP,nt,++nsnap,SNAP_FORMAT,SNAP,vx,vy,vz,sxx,syy,szz,u,pi,
					      IDX,IDY,IDZ,1,1,1,NX,NY,NZ);
				    lsnap=lsnap+iround(TSNAPINC/DT);
				}			 

			} /* end of loop over timesteps forward propagation*/

			if (SEISMO){
			  	saveseis(FP,sectionvx,sectionvy,sectionvz,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos,ishot,ns,0,iteration);
			}	
			
			/* output timing information (real times for update and exchange) */
			/*if (LOG)
			if (MYID==0) timing(time_v_update,time_s_update, time_s_exchange,time_v_exchange,time_timestep, ishot);*/
		
		
			if(METHOD==1){
				MPI_Barrier(MPI_COMM_WORLD);
				if(METHOD){
					/* exchange frequency domain wavefields */
					exchange_Fv(Ffvx,Ffvy,Ffvz,nf, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, 
											bufferfro_to_bac, bufferbac_to_fro, req_send, req_rec,0);
					exchange_Fv(Ffivx,Ffivy,Ffivz,nf, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, 
											bufferfro_to_bac, bufferbac_to_fro, req_send, req_rec,0);}
				
				/*---------------START INVERSION-------------*/
				

				if(ntr_loc>0){
					/* read seismic "observed" data (xcomp.) */
					readseis(ishot, sectionread, sectionreadf, ntr_loc, ns,1);
					if(FILT==1){
					filt_seis(sectionread,ntr_loc,NT,finv[nf-1]);}
					/* calculate residuals and L2 Norm */
					residual(sectionread, sectionreadf,sectionvx,sectionvxdiff,ntr_loc,ns,&L2,&L2f);
									
					/* read seismic "observed" data (ycomp.) */
					readseis(ishot, sectionread, sectionreadf, ntr_loc, ns,2);
					if(FILT==1){
					filt_seis(sectionread,ntr_loc,NT,finv[nf-1]);}
					/* calculate residuals and L2 Norm */
					residual(sectionread, sectionreadf,sectionvy,sectionvydiff,ntr_loc,ns,&L2,&L2f);
					
					/* read seismic "observed" data (zcomp.)*/
					readseis(ishot, sectionread, sectionreadf, ntr_loc, ns,3);
					if(FILT==1){
					filt_seis(sectionread,ntr_loc,NT,finv[nf-1]);}
					
					/* calculate residuals and L2 Norm */
					residual(sectionread, sectionreadf,sectionvz,sectionvzdiff,ntr_loc,ns,&L2,&L2f);
					
				}	
					
				MPI_Barrier(MPI_COMM_WORLD);
				buf=L2;
				MPI_Allreduce(&buf,&L2,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);	
				L2all+=L2;
						
				for (ishot1=1;ishot1<=NSHOTS_STEP;ishot1++){
				  /*printf("ishot1=%d",ishot1);*/
					if(ishot==((nshots/NSHOTS_STEP)*(ishot1-1)+1)){
						misfit[0]+=L2;}
				}
				if(MYID==0) fprintf(FI,"\n L2=%e, L2all=%e, misfit[0]=%e\n", L2, L2all,misfit[0]);
				
				
			/*ishot1=(nshots/NSHOTS_STEP)*(ishot-1)+1;*/
				
				/*---------------------------- backpropagation----------------------------------------------*/
				
				fprintf(FP,"\n\n *********** STARTING TIME STEPPING BACK PROPAGATION ***************\n");		  
				pshot_loc=1;
				for(pshot=0;pshot<=hloop;pshot++){
				
				  pshot1=pshot;
				lsamp=NDTSHIFT+1;
				nlsamp=1;
				if(pshot>0){/*only for Hessian wavefields*/
					pshot1=pshot;
					ntr_loc=rnum_loc[pshot];
					fprintf(FP,"\nstart Hessian wavefields pshot %d, ntr=%d\n", pshot,ntr);
					
					if(ntr_loc>0){
						for(i=1;i<=4;i++) {srcpos_loc_back[i][1]=recpos_loc[i][pshot_loc];srcpos_loc_back[6][1]=1.0;}
						wavelet(srcpos_loc_back,ntr_loc,5,sectionvydiff);
						if(FILT)filt_seis(sectionvydiff,ntr_loc,NT,finv[nf-1]);
						for(i=1;i<=NT;i++){
							sectionvxdiff[1][i]=0.0;
							sectionvzdiff[1][i]=0.0;
						}
						pshot_loc++;
					}
				}

				/* Initialisieren de Wellenfeldes mit Nullen */
				zero_wavefield(NX,NY,NZ,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rxx,ryy,rzz,rxy,ryz,rxz,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_sxy_y,psi_syy_y,psi_syz_y,psi_sxz_z,psi_syz_z,psi_szz_z,psi_vxx,psi_vyx,psi_vzx,psi_vxy,psi_vyy,psi_vzy,psi_vxz,psi_vyz,psi_vzz);
				
				if(MYID==0) fprintf(FP,"\n ****************************************\n ");
				
				for (nt=1;nt<=NT;nt++){
				  if(MYID==0) if(!(nt%(NT/40))) fprintf(FP,"*");

					time_v_update[nt]=0.0;
					time_s_update[nt]=0.0;

					/* Check if simulation is still stable */
					if (isnan(vy[NY/2][NX/2][NZ/2])) err(" Simulation is unstable !"); /* maybe just breaking the loop would be better */
		
					/* update of particle velocities */	
					time_v_update[nt]+=update_v(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],nt,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rho,rjp,rkp,rip,srcpos_loc_back,sectionvxdiff,sectionvydiff,sectionvzdiff,ntr_loc,absorb_coeff,1);

					if(ABS_TYPE==1){update_v_CPML(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],nt,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rho,rjp,rkp,rip,srcpos_loc_back,sectionvxdiff,sectionvydiff,sectionvzdiff,ntr_loc,absorb_coeff,1,K_x,a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,K_z,a_z,b_z,K_z_half,a_z_half,b_z_half,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_sxy_y,psi_syy_y,psi_syz_y,psi_sxz_z,psi_syz_z,psi_szz_z);}
			
					/* exchange values of particle velocities at grid boundaries between PEs */
					time_v_exchange[nt]+=exchange_v(vx,vy,vz, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, bufferfro_to_bac, bufferbac_to_fro, req_send, req_rec);
					/*time_s_update[nt]+=update_s(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rxx,ryy,rzz,rxy,ryz,rxz,pi,u,uipjp,ujpkp,uipkp,taus,tausipjp,tausjpkp,tausipkp,taup,eta);
					if(ABS_TYPE==1){
					update_s_CPML(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rxx,ryy,rzz,rxy,ryz,rxz,pi,u,uipjp,ujpkp,uipkp,taus,tausipjp,tausjpkp,tausipkp,taup,eta,K_x,a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,K_z,a_z,b_z,K_z_half,a_z_half,b_z_half,psi_vxx,psi_vyx,psi_vzx,psi_vxy,psi_vyy,psi_vzy,psi_vxz,psi_vyz,psi_vzz);}*/

					time_s_update[nt]+=update_s_elastic(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,pi,u,uipjp,ujpkp,uipkp);

					if(ABS_TYPE==1){				
					  update_s_CPML_elastic(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,pi,u,uipjp,ujpkp,uipkp,K_x,a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,K_z,a_z,b_z,K_z_half,a_z_half,b_z_half,psi_vxx,psi_vyx,psi_vzx,psi_vxy,psi_vyy,psi_vzy,psi_vxz,psi_vyz,psi_vzz);}
					
					time_s_exchange[nt]+=exchange_s(sxx,syy,szz,sxy,syz,sxz,sbufferlef_to_rig,sbufferrig_to_lef,sbuffertop_to_bot,sbufferbot_to_top,sbufferfro_to_bac,sbufferbac_to_fro, sreq_send, sreq_rec);
					  
					if(nt%ntast==0&&(METHOD)){
						discfourier(1,NX,1,NY,1,NZ,nt,vx,vy,vz,Fbvx,Fbvy,Fbvz, Fbivx, Fbivy,Fbivz,finv,nf,ntast,pshot1,1);
						l++;
					}
					
					/* stress free surface ? */
					if ((FREE_SURF) && (POS[2]==0))
						/*surface(1,u,pi,taus,taup,eta,sxx,syy,szz,sxy,syz,rxx,ryy,rzz,vx,vy,vz);*/		
						surface_elastic(1,u,pi,sxx,syy,szz,sxy,syz,vx,vy,vz);
					/* store amplitudes at receivers in sectionvx-sectionvz */
					if ((SEISMO) && (ntr>0) && (nt==lsamp)){
						seismo(nlsamp,ntr,recpos_loc,sectionvx,sectionvy,sectionvz,
							sectiondiv,sectioncurl,sectionp,vx,vy,vz,sxx,syy,szz,pi,u);
						nlsamp++;
						lsamp+=NDT;
					}    		

					/* save snapshot in file (backpropagation) */
					/*if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){
					    snap(FP,nt,++nsnap,SNAP_FORMAT,SNAP,vx,vy,vz,sxx,syy,szz,u,pi,
						      IDX,IDY,IDZ,1,1,1,NX,NY,NZ);
					    lsnap=lsnap+iround(TSNAPINC/DT);
					}	*/			

				}/*end of time loop*/ 

				fprintf(FP,"\n End backpropagation \n");
				/*saveseis(FP,sectionvx,sectionvy,sectionvz,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos,ishot,ns,0,pshot+10 );*/
				/*saveseis(FP,sectionvxdiff,sectionvydiff,sectionvzdiff,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos,ishot,ns,0,pshot+100 );*/
				MPI_Barrier(MPI_COMM_WORLD);
				
				if(!pshot){
					exchange_Fv(Fbvx,Fbvy,Fbvz,nf, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, 
											bufferfro_to_bac, bufferbac_to_fro, req_send, req_rec,0);
					exchange_Fv(Fbivx,Fbivy,Fbivz,nf, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, 
											bufferfro_to_bac, bufferbac_to_fro, req_send, req_rec,0);						 
					MPI_Barrier(MPI_COMM_WORLD);
								
					gradient_F(NX,NY,NZ,Ffvx,Ffvy,Ffvz,Ffivx,Ffivy,Ffivz,Fbvx,Fbvy,Fbvz,Fbivx,Fbivy,Fbivz,grad1,grad2,grad3,nt,rho,pi,u,finv,nf,iteration);
				}
				} /*hloop*/
				if(iteration==1 && HESS && !READ_HESS){
					exchange_Fv(Fbvx,Fbvy,Fbvz,nf, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, 
												bufferfro_to_bac, bufferbac_to_fro, req_send, req_rec,ntr_hess);
					exchange_Fv(Fbivx,Fbivy,Fbivz,nf, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, 
												bufferfro_to_bac, bufferbac_to_fro, req_send, req_rec,ntr_hess);						 
					MPI_Barrier(MPI_COMM_WORLD);
					hess_F(NX,NY,NZ,Ffvx,Ffvy,Ffvz,Ffivx,Ffivy,Ffivz,Fbvx,Fbvy,Fbvz, Fbivx,Fbivy,Fbivz, hess1,hess2,hess3,nt,rho, pi, u, finv, nf,ntr_hess);
				}
				hloop=0; ntr_loc=ntr;
				if(iteration==1 && HESS && !READ_HESS){
				for(i=1;i<=4;i++){
				for(j=1;j<=ntr_loc;j++){srcpos_loc_back[i][j]=recpos_loc[i][j];	}}}
				
			} /*end if(METHOD)*/
		}/*end of loop over shots */	
		MPI_Barrier(MPI_COMM_WORLD);
		
		if(METHOD){
			/*output Hessian*/
			if(HESS&&iteration==1) outgrad(NX,NY,NZ,hess1,hess2,hess3,finv[0],iteration, HESS_FILE);
		
			/*output "raw" gradient*/
			fprintf(FP,"\n raw Gradient: \n");
			outgrad(NX,NY,NZ,grad1,grad2,grad3,finv[0],iteration, GRAD_FILE);
			
			if(HESS) hess_apply(1,NX,1,NY,1,NZ,grad1,grad2,grad3,hess1,hess2,hess3,finv[0],iteration);
			
			/*preconditioning of gradient*/
			precon_grad(NX,NY,NZ,grad1,grad2,grad3,nsrc,srcpos,ntr_glob,recpos,finv[0],iteration,cdf);
			outgrad(NX,NY,NZ,grad1,grad2,grad3,finv[0],iteration+1000, GRAD_FILE);
			if(LBFGS){
				if(it_group>1){
					/*lbfgs(grad1, hess1, bfgsscale1, bfgsmod1, bfgsgrad1,iteration);*/
					lbfgs(grad1, grad2, grad3, bfgsscale1, bfgsmod1, bfgsgrad1,it_group);
					/*lbfgs(grad3, hess3, bfgsscale3, bfgsmod3, bfgsgrad3,iteration);*/
				}
				else lbfgs_savegrad(grad1,grad2,grad3,bfgsgrad1);
			}
			if(!LBFGS) conjugate(NX,NY,NZ,grad1,grad2,grad3,gradprior1,gradprior2,gradprior3,gradprior4,gradprior5,gradprior6,beta,iteration,cdf);
			outgrad(NX,NY,NZ,grad1,grad2,grad3,finv[0],iteration+2000, GRAD_FILE);
		
		/*---------------------------------------steplength calculation----------------------------------------------------------------*/
		
			fprintf(FP,"\n\n *********** STEPLENGTH CALCULATION ITERATION %d *********** \n", iteration);
				if(LBFGS||!LBFGS){
			for (steptest=1;steptest<=2;steptest++){
				if(cdf==1) step[4]=TESTSTEP;
				
				step[steptest]=0.0; step[steptest]=step[4]*steptest;
			
				fprintf(FP,"\n\n %d. test steplength: step[%d]=%.2e \n",steptest,steptest,step[steptest]);
				fprintf(FP," ------------------------------------ \n");
				
				if(MYID==0) fprintf(FI,"\n steptest %d: steplength=%e \n",steptest ,step[steptest]);
				
				cpmodel(NX,NY,NZ,rho,pi,u,testrho,testpi,testu);
				modelupdate(NX,NY,NZ,grad1,grad2,grad3,testrho,testpi,testu,bfgsmod1, step[steptest],beta,it_group);
				matcopy(testrho,testpi,testu,taus,taup);
				if (FREE_SURF)constant_boundary(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],testpi, testu, testrho);
				av_mat(testrho,testpi,testu,taus,taup,uipjp,ujpkp,uipkp,tausipjp,tausjpkp,tausipkp,rjp,rkp,rip);		
				
				L2=0.0;
				L2f=0.0;
			
				for (ishot=1;ishot<=NSHOTS_STEP;ishot++){
				  
					ishot1=(nshots/NSHOTS_STEP)*(ishot-1)+1;
				
					fprintf(FP,"\n **********  Starting simulation for shot %d of %d  ********** \n",ishot1,nshots);
								
					if (RUN_MULTIPLE_SHOTS){
						nsrc_loc=snum_loc[ishot1];
						if(nsrc_loc>0){
							for (nt=4;nt<=7;nt++) srcpos_loc[nt][1]=srcpos[nt][ishot1];
							srcpos_loc[1][1]=(float)(((iround(srcpos[1][ishot1]/DX)-1)%IENDX)+1);
							srcpos_loc[2][1]=(float)(((iround(srcpos[2][ishot1]/DY)-1)%IENDY)+1);
							srcpos_loc[3][1]=(float)(((iround(srcpos[3][ishot1]/DZ)-1)%IENDZ)+1);
							}
					}
							
					if(nsrc_loc>0){	
						wavelet(srcpos_loc,nsrc_loc,SOURCE_SHAPE,signals);	
						if(FILT==1) filt_seis(signals,nsrc_loc,NT,finv[nf-1]);}
					zero_wavefield(NX,NY,NZ,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rxx,ryy,rzz,rxy,ryz,rxz,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_sxy_y,psi_syy_y,psi_syz_y,psi_sxz_z,psi_syz_z,psi_szz_z,psi_vxx,psi_vyx,psi_vzx,psi_vxy,psi_vyy,psi_vzy,psi_vxz,psi_vyz,psi_vzz);
								
					lsamp=NDTSHIFT+1;
					nlsamp=1;
					
					
					if(MYID==0) fprintf(FP,"\n ****************************************\n ");
				
				
					
					for (nt=1;nt<=NT;nt++){
						if(MYID==0) if(!(nt%(NT/40))) fprintf(FP,"*");
						time_v_update[nt]=0.0;
						time_s_update[nt]=0.0;

						if (isnan(vy[NY/2][NX/2][NZ/2])) err(" Simulation is unstable !"); 
								
						time_v_update[nt]+=update_v(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],nt,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,testrho,rjp, rkp,rip,srcpos_loc,signals,signals,signals,nsrc_loc,absorb_coeff,0);

						if(ABS_TYPE==1){	update_v_CPML(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],nt,vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,testrho,rjp,rkp,rip,srcpos_loc,signals,signals,signals,nsrc_loc,absorb_coeff,0,K_x,a_x,b_x,K_x_half,a_x_half,b_x_half, K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,K_z,a_z,b_z,K_z_half,a_z_half,b_z_half,psi_sxx_x,psi_sxy_x,psi_sxz_x,psi_sxy_y,psi_syy_y,psi_syz_y,psi_sxz_z,psi_syz_z,psi_szz_z);}
		;
						time_v_exchange[nt]+=exchange_v(vx,vy,vz, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, bufferfro_to_bac, bufferbac_to_fro, req_send, req_rec);
			    
						/*time_s_update[nt]+=update_s(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rxx,ryy,rzz,rxy,ryz,rxz,testpi,testu,uipjp,ujpkp,uipkp,taus,tausipjp,tausjpkp,tausipkp,taup,eta);
						if(ABS_TYPE==1){
						update_s_CPML(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,rxx,ryy,rzz,rxy,ryz,rxz,testpi,testu,uipjp,ujpkp,uipkp,taus,tausipjp,tausjpkp,tausipkp,taup,eta,K_x,a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,K_z,a_z,b_z,K_z_half,a_z_half,b_z_half,psi_vxx,psi_vyx,psi_vzx,psi_vxy,psi_vyy,psi_vzy,psi_vxz,psi_vyz,psi_vzz);}*/

						time_s_update[nt]+=update_s_elastic(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,testpi,testu,uipjp,ujpkp,uipkp);

						if(ABS_TYPE==1){	  update_s_CPML_elastic(xb[0],xb[1],yb[0],yb[1],zb[0],zb[1],vx,vy,vz,sxx,syy,szz,sxy,syz,sxz,testpi,testu,uipjp,ujpkp,uipkp,K_x,a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,K_z,a_z,b_z,K_z_half,a_z_half,b_z_half,psi_vxx,psi_vyx,psi_vzx,psi_vxy,psi_vyy,psi_vzy,psi_vxz,psi_vyz,psi_vzz);}
		
						psource(nt,sxx,syy,szz,srcpos_loc,signals,nsrc_loc);
						
						time_s_exchange[nt]+=exchange_s(sxx,syy,szz,sxy,syz,sxz,sbufferlef_to_rig,sbufferrig_to_lef,sbuffertop_to_bot,sbufferbot_to_top,sbufferfro_to_bac,sbufferbac_to_fro, sreq_send, sreq_rec);

						if ((FREE_SURF) && (POS[2]==0))
							/*surface(1,testu,testpi,taus,taup,eta,sxx,syy,szz,sxy,syz,rxx,ryy,rzz,vx,vy,vz);*/		
							surface_elastic(1,u,pi,sxx,syy,szz,sxy,syz,vx,vy,vz);
						if ((SEISMO) && (ntr>0) && (nt==lsamp)){
							seismo(nlsamp,ntr,recpos_loc,sectionvx,sectionvy,sectionvz,sectiondiv,sectioncurl,sectionp,vx,vy,vz,sxx,syy,szz,testpi,testu);
							nlsamp++;
							lsamp+=NDT;
						}    		
					} 
						
					if(ntr_loc>0){
						readseis(ishot1, sectionread, sectionreadf, ntr_loc, ns,1);
						if(FILT==1){filt_seis(sectionread,ntr_loc,NT,finv[nf-1]);}
						residual(sectionread, sectionreadf,sectionvx,sectionvxdiff,ntr_loc,ns,&L2,&L2f); 
						readseis(ishot1, sectionread, sectionreadf, ntr_loc, ns,2);
						if(FILT==1){filt_seis(sectionread,ntr_loc,NT,finv[nf-1]);}
						residual(sectionread, sectionreadf,sectionvy,sectionvydiff,ntr_loc,ns,&L2,&L2f); 
						readseis(ishot1, sectionread, sectionreadf, ntr_loc, ns,3);
						if(FILT==1){filt_seis(sectionread,ntr_loc,NT,finv[nf-1]);}
						residual(sectionread, sectionreadf,sectionvz,sectionvzdiff,ntr_loc,ns,&L2,&L2f); 	
						}
		
				} /*ishot*/
				
				MPI_Barrier(MPI_COMM_WORLD);
				buf=L2;
				MPI_Allreduce(&buf,&L2,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);	
				
				misfit[steptest]=0.0;
				misfit[steptest]=L2;
				if(MYID==0)fprintf(FI,"\n L2=%e",misfit[steptest]);
				
			} /*steptest*/
			

			if(MYID==0)fprintf(FP,"\n\n Steplength Parabel \n");
			MPI_Barrier(MPI_COMM_WORLD);
			steplength(misfit,step,iteration, it_group); /*find optimal steplength*/
			}
			if(MYID==0)fprintf(FP,"\n stepength calculation finished\n");
			
			
			modelupdate(NX,NY,NZ,grad1,grad2,grad3,rho,pi,u,bfgsmod1,step[3],beta,it_group);
					
			if(MYID==0)fprintf(FP,"\n Modeloutput \n");
			outmod(NX,NY,NZ,rho,pi,u,iteration);
			
			if(MYID==0){
				time4=MPI_Wtime();
				fprintf(FP,"\n Iteration finished;%4.2f\n", time4-time2);
			}
		}/*fwi loop*/
	}/*iteration loop*/
	
	fprintf(FP,"\n Inversion finished \n");
	if(MYID==0)fprintf(FI,"\n\n *********** INVERSION FINISHED*********** \n\n");

	/*de-allocation of memory */
	l=1;
	if(ABS_TYPE==1 && FDORDER==2){l=2;}
		
	if(POS[2]==0){
	free_f3tensor(vx,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	free_f3tensor(vy,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	free_f3tensor(vz,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
        free_f3tensor(sxy,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	free_f3tensor(syz,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	if(METHOD){
		free_f4tensor(Ffvx,1,NFMAX,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Ffvy,1,NFMAX,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Ffvz,1,NFMAX,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Fbvx,1,NFMAX*(ntr_hess+1),0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Fbvy,1,NFMAX*(ntr_hess+1),0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Fbvz,1,NFMAX*(ntr_hess+1),0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Ffivx,1,NFMAX,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Ffivy,1,NFMAX,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Ffivz,1,NFMAX,0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Fbivx,1,NFMAX*(ntr_hess+1),0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Fbivy,1,NFMAX*(ntr_hess+1),0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Fbivz,1,NFMAX*(ntr_hess+1),0-FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	}
	}	

	if(POS[2]>0){
	free_f3tensor(vx,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	free_f3tensor(vy,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	free_f3tensor(vz,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
        free_f3tensor(sxy,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	free_f3tensor(syz,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	/*free_f4tensor(fvx,1,NT,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	free_f4tensor(fvy,1,NT,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	free_f4tensor(fvz,1,NT,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);*/
	if(METHOD){
		free_f4tensor(Ffvx,1,NFMAX,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Ffvy,1,NFMAX,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Ffvz,1,NFMAX,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Fbvx,1,NFMAX*(ntr_hess+1),1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Fbvy,1,NFMAX*(ntr_hess+1),1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Fbvz,1,NFMAX*(ntr_hess+1),1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Ffivx,1,NFMAX,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Ffivy,1,NFMAX,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Ffivz,1,NFMAX,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Fbivx,1,NFMAX*(ntr_hess+1),1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Fbivy,1,NFMAX*(ntr_hess+1),1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
		free_f4tensor(Fbivz,1,NFMAX*(ntr_hess+1),1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	}
	}	

	
	free_f3tensor(sxz,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	free_f3tensor(sxx,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	free_f3tensor(syy,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
	free_f3tensor(szz,1-l*FDORDER/2,NY+l*FDORDER/2,1-l*FDORDER/2,NX+l*FDORDER/2,1-l*FDORDER/2,NZ+l*FDORDER/2);
        
	if(METHOD){
		free_f3tensor(grad1,1,NX,1,NY,1,NZ);
		free_f3tensor(grad2,1,NX,1,NY,1,NZ);
		free_f3tensor(grad3,1,NX,1,NY,1,NZ);
		free_f3tensor(gradprior1,1,NX,1,NY,1,NZ);
		free_f3tensor(gradprior2,1,NX,1,NY,1,NZ);
		free_f3tensor(gradprior3,1,NX,1,NY,1,NZ);
		free_f3tensor(gradprior4,1,NX,1,NY,1,NZ);
		free_f3tensor(gradprior5,1,NX,1,NY,1,NZ);
		free_f3tensor(gradprior6,1,NX,1,NY,1,NZ);
		if(HESS){
			free_f3tensor(hess1,1,NY,1,NX,1,NZ);
			free_f3tensor(hess2,1,NY,1,NX,1,NZ);
			free_f3tensor(hess3,1,NY,1,NX,1,NZ);
		}
	}
	
	if(LBFGS){
		free_matrix(bfgsmod1,1,BFGSNUM,1,NUMPAR*NX*NY*NZ);
		free_matrix(bfgsgrad1,1,BFGSNUM,1,NUMPAR*NX*NY*NZ);
		free_vector(bfgsscale1,1,BFGSNUM);
	}
	
	if(ABS_TYPE==1){
	  
	free_vector(K_x,1,2*FW);
	free_vector(alpha_prime_x,1,2*FW);
	free_vector(a_x,1,2*FW);
	free_vector(b_x,1,2*FW);
	free_vector(K_x_half,1,2*FW);
	free_vector(alpha_prime_x_half,1,2*FW);
	free_vector(a_x_half,1,2*FW);
	free_vector(b_x_half,1,2*FW);

	free_vector(K_y,1,2*FW);
	free_vector(alpha_prime_y,1,2*FW);
	free_vector(a_y,1,2*FW);
	free_vector(b_y,1,2*FW);
	free_vector(K_y_half,1,2*FW);
	free_vector(alpha_prime_y_half,1,2*FW);
	free_vector(a_y_half,1,2*FW);
	free_vector(b_y_half,1,2*FW);

	free_vector(K_z,1,2*FW);
	free_vector(alpha_prime_z,1,2*FW);
	free_vector(a_z,1,2*FW);
	free_vector(b_z,1,2*FW);
	free_vector(K_z_half,1,2*FW);
	free_vector(alpha_prime_z_half,1,2*FW);
	free_vector(a_z_half,1,2*FW);
	free_vector(b_z_half,1,2*FW);

	free_f3tensor(psi_sxx_x,1,NY,1,2*FW,1,NZ);
	free_f3tensor(psi_syy_y,1,2*FW,1,NX,1,NZ);
	free_f3tensor(psi_szz_z,1,NY,1,NX,1,2*FW);
	free_f3tensor(psi_sxy_x,1,NY,1,2*FW,1,NZ);
	free_f3tensor(psi_sxy_y,1,2*FW,1,NX,1,NZ);
	free_f3tensor(psi_sxz_x,1,NY,1,2*FW,1,NZ);
	free_f3tensor(psi_sxz_z,1,NY,1,NX,1,2*FW);
	free_f3tensor(psi_syz_y,1,2*FW,1,NX,1,NZ);
	free_f3tensor(psi_syz_z,1,NY,1,NX,1,2*FW);
	
	free_f3tensor(psi_vxx,1,NY,1,2*FW,1,NZ);
	free_f3tensor(psi_vyy,1,2*FW,1,NX,1,NZ);
	free_f3tensor(psi_vzz,1,NY,1,NX,1,2*FW);
	free_f3tensor(psi_vxy,1,2*FW,1,NX,1,NZ);
	free_f3tensor(psi_vyx,1,NY,1,2*FW,1,NZ);
	free_f3tensor(psi_vxz,1,NY,1,NX,1,2*FW);
	free_f3tensor(psi_vzx,1,NY,1,2*FW,1,NZ);
 	free_f3tensor(psi_vyz,1,NY,1,NX,1,2*FW);
	free_f3tensor(psi_vzy,1,2*FW,1,NX,1,NZ);
	}

	if(L){
	free_f3tensor(rxx,1,NY,1,NX,1,NZ);
	free_f3tensor(ryy,1,NY,1,NX,1,NZ);
	free_f3tensor(rzz,1,NY,1,NX,1,NZ);
	free_f3tensor(rxy,1,NY,1,NX,1,NZ);
	free_f3tensor(ryz,1,NY,1,NX,1,NZ);
	free_f3tensor(rxz,1,NY,1,NX,1,NZ);
	free_f3tensor(taus,0,NY+1,0,NX+1,0,NZ+1);
	free_f3tensor(taup,0,NY+1,0,NX+1,0,NZ+1);
	free_vector(eta,1,L);
	free_f3tensor(tausipjp,1,NY,1,NX,1,NZ);
	free_f3tensor(tausjpkp,1,NY,1,NX,1,NZ);
	free_f3tensor(tausipkp,1,NY,1,NX,1,NZ);
	}

	free_f3tensor(rho,0,NY+1,0,NX+1,0,NZ+1);
	free_f3tensor(pi,0,NY+1,0,NX+1,0,NZ+1);
	free_f3tensor(u,0,NY+1,0,NX+1,0,NZ+1);
	free_f3tensor(absorb_coeff,1,NY,1,NX,1,NZ);
	free_f3tensor(testrho ,0,NY+1,0,NX+1,0,NZ+1);
	free_f3tensor(testpi,0,NY+1,0,NX+1,0,NZ+1);
	free_f3tensor(testu,0,NY+1,0,NX+1,0,NZ+1);

	/* averaged material parameters */
	free_f3tensor(uipjp,1,NY,1,NX,1,NZ);
	free_f3tensor(ujpkp,1,NY,1,NX,1,NZ);
	free_f3tensor(uipkp,1,NY,1,NX,1,NZ);
	free_f3tensor(rjp,1,NY,1,NX,1,NZ);
	free_f3tensor(rkp,1,NY,1,NX,1,NZ);
	free_f3tensor(rip,1,NY,1,NX,1,NZ);


	
	free_f3tensor(bufferlef_to_rig,1,NY,1,NZ,1,nf1);
	free_f3tensor(bufferrig_to_lef,1,NY,1,NZ,1,nf2);
	free_f3tensor(buffertop_to_bot,1,NX,1,NZ,1,nf1);
	free_f3tensor(bufferbot_to_top,1,NX,1,NZ,1,nf2);
	free_f3tensor(bufferfro_to_bac,1,NY,1,NX,1,nf1);
	free_f3tensor(bufferbac_to_fro,1,NY,1,NX,1,nf2);

	free_f3tensor(sbufferlef_to_rig,1,NY,1,NZ,1,nf2);
	free_f3tensor(sbufferrig_to_lef,1,NY,1,NZ,1,nf1);
	free_f3tensor(sbuffertop_to_bot,1,NX,1,NZ,1,nf2);
	free_f3tensor(sbufferbot_to_top,1,NX,1,NZ,1,nf1);
	free_f3tensor(sbufferfro_to_bac,1,NY,1,NX,1,nf2);
	free_f3tensor(sbufferbac_to_fro,1,NY,1,NX,1,nf1);

	/* free memory for global source positions */
	free_imatrix(recpos,1,3,1,ntr_glob);

	/* free memory for source positions */
	if (nsrc_loc>0){	
		if(RUN_MULTIPLE_SHOTS){	free_matrix(signals,1,1,1,NT);
					free_matrix(srcpos_loc,1,7,1,1);}
		else {			free_matrix(signals,1,nsrc_loc,1,NT);
					free_matrix(srcpos_loc,1,7,1,nsrc_loc);}
	}
        if(METHOD)free_matrix(srcpos_loc_back,1,7,1,ntr);
	free_matrix(srcpos,1,7,1,nsrc);
	free_matrix(srcpos1,1,6,1,1);
	free_ivector(snum_loc,1,nsrc);
	free_ivector(rnum_loc,1,ntr);

	if ((ntr>0) && (SEISMO)){	

      		free_imatrix(recpos_loc,1,4,1,ntr);
		switch (SEISMO){
		case 1 : /* particle velocities only */
			free_matrix(sectionvx,1,ntr,1,ns);
			free_matrix(sectionvy,1,ntr,1,ns);		
			free_matrix(sectionvz,1,ntr,1,ns);		
			break;	
		case 2 : /* pressure only */
			free_matrix(sectionp,1,ntr,1,ns);
			break;	
		case 3 : /* curl and div only */
			free_matrix(sectioncurl,1,ntr,1,ns);
			free_matrix(sectiondiv,1,ntr,1,ns);
			break;	
		case 4 : /* everything */
			free_matrix(sectionvx,1,ntr,1,ns);
			free_matrix(sectionvy,1,ntr,1,ns);
			free_matrix(sectionvz,1,ntr,1,ns);
			free_matrix(sectionp,1,ntr,1,ns);
			free_matrix(sectioncurl,1,ntr,1,ns);
			free_matrix(sectiondiv,1,ntr,1,ns);		
			break;
		}	

	}	

	
	/* free inversion variables */
	free_matrix(sectionread,1,ntr_glob,1,ns);
	free_matrix(sectionreadf,1,ntr_glob,1,1);
	free_matrix(sectionvxdiff,1,ntr_glob,1,ns);
	free_matrix(sectionvydiff,1,ntr_glob,1,ns);
	free_matrix(sectionvzdiff,1,ntr_glob,1,ns);
	
	/* de-allocate buffer for messages */
	MPI_Buffer_detach(buff_addr,&buffsize);

	MPI_Barrier(MPI_COMM_WORLD);

	/* merge snapshot files created by the PEs into one file */
	/* if ((SNAP) && (MYID==0)) snapmerge(nsnap);*/
	
   	free_ivector(xb,0,1);
	free_ivector(yb,0,1);
	free_ivector(zb,0,1);

	free_vector(misfit,0,3);
	free_vector(step,0,4);
	free_vector(beta,0,2);
	free_vector(finv,0,NFMAX-1);
	
        /* free timing arrays */
	free_dvector(time_v_update,1,NT);
	free_dvector(time_s_update,1,NT);
	free_dvector(time_s_exchange,1,NT);
	free_dvector(time_v_exchange,1,NT);
	free_dvector(time_timestep,1,NT);

	
	if (MYID==0){
		fprintf(FP,"\n **Info from main (written by PE %d): \n",MYID);
		time4=MPI_Wtime();
	        fprintf(FP," Total real time of program: %4.2f seconds.\n\n",time4-time1);
		fprintf(FP," ***********************************************************\n");
		fprintf(FP," IFOS3D has finished.\n");
		fprintf(FP," ***********************************************************\n\n");
	}
	

	fclose(FP);
	fclose(FI);

	MPI_Finalize();

	return 0;

}
