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

/* ---------------------------------------------------------------------
 *  Computation of local source coordinates (within each subgrid)
 -----------------------------------------------------------------------*/
 
 #include "fd.h"
float **splitsrc(float **srcpos,int *nsrc_loc, int nsrc, int *snum_loc)
{

	extern int IENDX, IENDY, IENDZ, MYID, POS[4], RUN_MULTIPLE_SHOTS,VERBOSE;
	extern float DX, DY, DZ;
	extern FILE *FP;

	int a,b,c,i=0,j,k;
	float ** srcpos_dummy, **srcpos_local=NULL;

	srcpos_dummy = fmatrix(1,7,1,nsrc);

	for (j=1;j<=nsrc;j++) {
		a=(iround(srcpos[1][j]/DX)-1)/IENDX;
		b=(iround(srcpos[2][j]/DY)-1)/IENDY;
		c=(iround(srcpos[3][j]/DZ)-1)/IENDZ;

		if ((POS[1]==a)&&(POS[2]==b)&&(POS[3]==c)) {
			i++;
			srcpos_dummy[1][i] = (float)(((iround(srcpos[1][j]/DX)-1)%IENDX)+1);
			srcpos_dummy[2][i] = (float)(((iround(srcpos[2][j]/DY)-1)%IENDY)+1);
			srcpos_dummy[3][i] = (float)(((iround(srcpos[3][j]/DZ)-1)%IENDZ)+1);
			srcpos_dummy[4][i] = srcpos[4][j];
			srcpos_dummy[5][i] = srcpos[5][j];
			srcpos_dummy[6][i] = srcpos[6][j];
			srcpos_dummy[7][i] = srcpos[7][j];; /* stype_loc might be scaled for minimum size by realloc later */
			snum_loc[j]        =1;
		}
		else 	snum_loc[j]        =0;
		
	}
	if(RUN_MULTIPLE_SHOTS) i=1;
	if (i>0) srcpos_local = fmatrix(1,7,1,i);
	for (k=1;k<=i;k++){
		srcpos_local[1][k] = srcpos_dummy[1][k];
		srcpos_local[2][k] = srcpos_dummy[2][k];
		srcpos_local[3][k] = srcpos_dummy[3][k];
		srcpos_local[4][k] = srcpos_dummy[4][k];
		srcpos_local[5][k] = srcpos_dummy[5][k];
		srcpos_local[6][k] = srcpos_dummy[6][k];
		srcpos_local[7][k] = srcpos_dummy[7][k];
	}

	free_matrix(srcpos_dummy,1,7,1,nsrc);
	if (!RUN_MULTIPLE_SHOTS) {
	if(VERBOSE){
	fprintf(FP,"\n **Message from splitsrc:\n");
	fprintf(FP," Splitting of source positions from global to local grids finished.\n");
	fprintf(FP," MYID= %d \t \t no. of sources= %d\n",MYID,i);
        if(MYID==0){
	printf("\n **Message from splitsrc:\n");
	printf(" Table of local source positions (in gridpoints), time-shift, centre frequency, amplitude and source type:\n");
	printf(" MYID\t  x\t  y\t  z\t  tshift  fc\t  amp\t stype\n");}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	for (j=1;j<=i;j++)
		printf(" %3d\t%4.0f\t%4.0f\t%4.0f\t%6.2f\t%6.2f\t%6.2f\t%4.2f\n",
		    MYID,srcpos_local[1][j],srcpos_local[2][j],srcpos_local[3][j],
				   srcpos_local[4][j],srcpos_local[5][j],srcpos_local[6][j],srcpos_local[7][j]);}
        MPI_Barrier(MPI_COMM_WORLD);

        *nsrc_loc=i;
	return srcpos_local;

}
