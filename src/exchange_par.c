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
 *   Exchange FD-Parameters between PEs                         
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void exchange_par(void){
	/* declaration of extern variables */
	extern int   NX, NY, NZ, SOURCE_SHAPE, SOURCE_TYPE, SNAP, SNAP_FORMAT, SNAP_PLANE;
	extern int DRX, DRZ, L, SRCREC, FDORDER;
	extern float DX, DY, DZ, TIME, DT, *FL, TS, TAU, PLANE_WAVE_DEPTH, PHI;
	extern float XREC1, XREC2, YREC1, YREC2, ZREC1, ZREC2;
	extern float ALPHA, BETA, VPPML;
	extern float REC_ARRAY_DEPTH, REC_ARRAY_DIST;
	extern int SEISMO, NDT, NDTSHIFT, NGEOPH, SEIS_FORMAT, FREE_SURF, READMOD, READREC;
	extern int BOUNDARY, REC_ARRAY, LOG, IDX, IDY, IDZ, ABS_TYPE;
	extern float TSNAP1, TSNAP2, TSNAPINC, FW, REFREC[4], DAMPING, FPML;
	extern char  MFILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE];
	extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE], SEIS_FILE[STRING_SIZE];
	extern char GRAD_FILE[STRING_SIZE],SEIS_OBS_FILE[STRING_SIZE],INV_FILE[STRING_SIZE];
	extern char MOD_OUT_FILE[STRING_SIZE], HESS_FILE[STRING_SIZE];
	extern int NPROC,NPROCX,NPROCY,NPROCZ, MYID, RUN_MULTIPLE_SHOTS, FDCOEFF;
	extern char  FILEINP[STRING_SIZE];
	extern int   LITTLEBIG, ASCIIEBCDIC, IEEEIBM;
	extern int METHOD;
	extern int ITMIN, ITMAX, FILT, NFMAX, TAST, NSHOTS_STEP, DAMPTYPE, HESS, READ_HESS, REC_HESS, EXTOBS;
	extern int BFGSNUM, NUMPAR, LBFGS;
	/*extern float F_INV;*/
	extern float TESTSTEP, WATER_HESS[3], WEIGHT[3], VP0, VS0, RHO0;
	int idum[57];
	float fdum[42];
	
	
	if (MYID == 0){ 
		fdum[1]  = DX;
                fdum[2]  = DY;
                fdum[3]  = DZ;
		fdum[4]  = TIME;
		fdum[5]  = DT;
		fdum[6]  = TS;
		fdum[7]  = PHI;
		fdum[8]  = 0.0;
		fdum[9]  = 0.0;
		fdum[10]  = TAU;
		fdum[11]  = FW;
		fdum[12]  = TSNAP1;
		fdum[13]  = TSNAP2;
		fdum[14]  = TSNAPINC;
		fdum[15]  = REFREC[1];
		fdum[16]  = REFREC[2];
		fdum[17]  = REFREC[3];
		fdum[18]  = XREC1;
		fdum[19]  = YREC1;
		fdum[20]  = ZREC1;
		fdum[21]  = XREC2;
		fdum[22]  = YREC2;
		fdum[23]  = ZREC2;
		fdum[24]  = DAMPING;
		fdum[25]  = FPML;
		fdum[26]  = REC_ARRAY_DEPTH;
		fdum[27]  = REC_ARRAY_DIST;
		fdum[28]  = PLANE_WAVE_DEPTH;
		fdum[29]  = ALPHA;
		fdum[30]  = BETA;
		fdum[31]  = VPPML;
		fdum[32]  = VP0;
		fdum[33]  = VS0;
		fdum[34]  = RHO0;
		fdum[35]  = WEIGHT[0];
		fdum[36]  = WEIGHT[1];
		fdum[37]  = WEIGHT[2];
		fdum[38]  = TESTSTEP;
		fdum[39]  = WATER_HESS[0];
		fdum[40]  = WATER_HESS[1];
		fdum[41]  = WATER_HESS[2];

		idum[0]  = FDORDER;
		idum[1]  = NPROCX;
		idum[2]  = NPROCY;
		idum[3]  = NPROCZ;
		idum[4]  = NPROCX*NPROCY*NPROCZ;
		idum[5]  = NX;
		idum[6]  = NY;
		idum[7]  = NZ;
		idum[8]  = SOURCE_SHAPE;
		idum[9]  = SOURCE_TYPE;
		idum[10]  = READMOD;
		idum[11]  = L;
		idum[12]  = FREE_SURF;
		idum[13]  = SNAP;
		idum[14]  = DRX;
		idum[15]  = DRZ;
		idum[16]  = BOUNDARY;
		idum[17]  = REC_ARRAY;
		idum[18]  = SRCREC;
		idum[19]  = LOG;
		idum[20]  = IDX;
		idum[21]  = IDY;
		idum[22]  = IDZ;
		idum[23]  = SNAP_FORMAT;
		idum[24]  = SEISMO;
		idum[25]  = READREC;
		idum[26]  = NGEOPH;
		idum[27]  = NDT;
		idum[28]  = NDTSHIFT;
		
		idum[29]  = SEIS_FORMAT;
		
		                                                                       
        	idum[35]  = RUN_MULTIPLE_SHOTS;                                                                       
		idum[36]  = SNAP_PLANE;
		idum[37]  = ABS_TYPE;
                idum[38]  = FDCOEFF;
		
		idum[39]  = ASCIIEBCDIC;
		idum[40]  = LITTLEBIG;
		idum[41]  = IEEEIBM;

		idum[42]  = METHOD;
		
		idum[43] = EXTOBS;
		idum[44] = ITMIN;
		idum[45] = ITMAX;
		idum[46] = FILT;
		idum[47] = NFMAX;
		idum[48] = TAST;
		idum[49] = NSHOTS_STEP;
		idum[50] = DAMPTYPE;
		idum[51] = HESS;
		idum[52] = READ_HESS;
		idum[53] = REC_HESS;;
		idum[54] = LBFGS;
		idum[55] = NUMPAR;
		idum[56] = BFGSNUM;
		}

	if (MYID != 0) FL=vector(1,L);
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(&idum,57,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&fdum,42,MPI_FLOAT,0,MPI_COMM_WORLD);

	MPI_Bcast(&SOURCE_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SIGNAL_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&MFILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SNAP_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&REC_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SEIS_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	
	MPI_Bcast(&FILEINP,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&GRAD_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&SEIS_OBS_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&INV_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&MOD_OUT_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&HESS_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
	

	MPI_Barrier(MPI_COMM_WORLD);

	DX=fdum[1];
	DY=fdum[2];
	DZ=fdum[3];
	TIME=fdum[4];
	DT=fdum[5];
	TS=fdum[6];
	PHI=fdum[7];

	TAU=fdum[10];
	FW=fdum[11];
	TSNAP1=fdum[12];
	TSNAP2=fdum[13];
	TSNAPINC=fdum[14];
	REFREC[1]=fdum[15];
	REFREC[2]=fdum[16];
	REFREC[3]=fdum[17];
	XREC1=fdum[18];
	YREC1=fdum[19];
	ZREC1=fdum[20];
	XREC2=fdum[21];
	YREC2=fdum[22];
	ZREC2=fdum[23];
	DAMPING=fdum[24];
	FPML=fdum[25];
	REC_ARRAY_DEPTH=fdum[26];
	REC_ARRAY_DIST=fdum[27];
	PLANE_WAVE_DEPTH=fdum[28];
	ALPHA = fdum[29];
	BETA = fdum[30];
        VPPML = fdum[31];
	VP0 = fdum[32];
	VS0 = fdum[33];
	RHO0 = fdum[34];
	WEIGHT[0] = fdum[35];
	WEIGHT[1] = fdum[36];
	WEIGHT[2] = fdum[37];
	TESTSTEP=fdum[38];
	WATER_HESS[0]=fdum[39];
	WATER_HESS[1]=	fdum[40];
	WATER_HESS[2]=	fdum[41];

	FDORDER = idum[0];
	NPROCX = idum[1];
	NPROCY = idum[2];
	NPROCZ = idum[3];
	NPROC  = idum[4];
	NX = idum[5];
	NY = idum[6];
	NZ = idum[7];
	SOURCE_SHAPE = idum[8];
	SOURCE_TYPE = idum[9];
	READMOD = idum[10];
	L = idum[11];
	FREE_SURF = idum[12];
	SNAP = idum[13];
	DRX = idum[14];
	DRZ = idum[15];
	BOUNDARY = idum[16];
	REC_ARRAY = idum[17];
	SRCREC = idum[18];
	LOG = idum[19];
	IDX = idum[20];
	IDY = idum[21];
	IDZ = idum[22];
	
	SNAP_FORMAT = idum[23];
	SEISMO = idum[24];
	READREC = idum[25];
	NGEOPH = idum[26];
	NDT = idum[27];
	NDTSHIFT = idum[28];
	
	SEIS_FORMAT = idum[29];
	
	RUN_MULTIPLE_SHOTS= idum[35];
	SNAP_PLANE= idum[36];
	ABS_TYPE= idum[37];
	FDCOEFF= idum[38];
	
	ASCIIEBCDIC=idum[39];
	LITTLEBIG=idum[40];
	IEEEIBM=idum[41];

	METHOD=idum[42];
	
	EXTOBS = idum[43];
	ITMIN = idum[44];
	ITMAX = idum[45];
	FILT = idum[46];
	NFMAX = idum[47];
	TAST = idum[48];
	NSHOTS_STEP = idum[49];
	DAMPTYPE = idum[50];
	HESS = idum[51];
	READ_HESS = idum[52];
	REC_HESS  =idum[53];
	LBFGS = idum[54];
	NUMPAR = idum[55];
	BFGSNUM = idum[56];
	
	MPI_Bcast(&FL[1],L,MPI_FLOAT,0,MPI_COMM_WORLD);

}

