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
 * definition of global variables used in IFOS3D
 * For the names of the global variables
 *  uppercase letters are used
 --------------------------------------------------------------------------*/


float DX, DY, DZ, TIME, DT, TS, PLANE_WAVE_DEPTH, PHI;
float TSNAP1, TSNAP2, TSNAPINC, *FL, TAU, REC_ARRAY_DEPTH, REC_ARRAY_DIST;
float XREC1, XREC2, YREC1, YREC2, ZREC1=0.0, ZREC2=0.0;
float REFREC[4]={0.0, 0.0, 0.0, 0.0}, DAMPING=8.0, VPPML, FPML;
int   SEISMO, NDT, NDTSHIFT, NGEOPH, SEIS_FORMAT, FREE_SURF, READMOD, MOD_FORMAT, READREC, REC_ARRAY, LOG, FDORDER, FW=0, ABS_TYPE, BLOCK;
int   NX, NY, NZ=1, NT, SOURCE_SHAPE, SOURCE_TYPE, SNAP, SNAP_FORMAT, BOUNDARY, SRCREC, SNAP_PLANE;
float ALPHA, BETA;
int   NXG, NYG, NZG, IDX, IDY, IDZ, L=1, NX1, NX2, NY1, NY2, NZ1, NZ2, DRX, DRZ, RUN_MULTIPLE_SHOTS, FDCOEFF;
char  SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE], INV_FILE[STRING_SIZE];
char  MFILE[STRING_SIZE], REC_FILE[STRING_SIZE];
char  SEIS_FILE[STRING_SIZE],MOD_FILE[STRING_SIZE];
int VERBOSE;

FILE  *FP=NULL;
FILE  *FI=NULL;
char  FILEINP[STRING_SIZE]; /* input file name (appears in SEG-Y header) */
int   LITTLEBIG, ASCIIEBCDIC, IEEEIBM; /* computer charcteristics */
int   FDMPIVERS; /* version of fdmpi 33: current 3D isotropic elastic (SSG), 32: current 3D isotropic acoustic (SSG) */

/* Mpi-variables */
int   NP, NPSP, NPROC, NPROCX, NPROCY, NPROCZ, MYID, IENDX, IENDY, IENDZ;
int   POS[4], INDEX[7];     
const int TAG1=1,TAG2=2, TAG3=3, TAG4=4, TAG5=5,TAG6=6; 


float FC,AMP, REFSRC[3], SRC_DT, SRCTSHIFT;
int SRC_MF, SIGNAL_FORMAT, SRCOUT_PAR[6], FSRC, JSRC, LSRC;
char SRCOUT_FILE[STRING_SIZE];

int METHOD;
float F_INV=0.0, TESTSTEP=0.0;
char MOD_OUT_FILE[STRING_SIZE], HESS_FILE[STRING_SIZE], GRAD_FILE[STRING_SIZE], SEIS_OBS_FILE[STRING_SIZE];
int DAMPTYPE, ITMIN, ITMAX, FILT, NFMAX, NSHOTS_STEP, TAST,EXTOBS;
int HESS, READ_HESS, REC_HESS, LBFGS, NUMPAR, BFGSNUM;
float WATER_HESS[3]={0.0, 0.0, 0.0}, WEIGHT[3]={0.0, 0.0, 0.0};
float VP0, VS0, RHO0;