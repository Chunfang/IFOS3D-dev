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

/*-----------------------------------------------
 * Find global grid positions for the receivers.
 -----------------------------------------------*/

#include "fd.h"

int **receiver(FILE *fp, int *ntr) {

	/* declaration of extern variables */
	extern char REC_FILE[STRING_SIZE];
	extern float XREC1, YREC1, ZREC1, XREC2, YREC2, ZREC2;
	extern float DX, DY, DZ, REFREC[4], REC_ARRAY_DEPTH, REC_ARRAY_DIST;
	extern int READREC, NGEOPH, NXG, NZG, REC_ARRAY, BOUNDARY;
	extern int MYID, DRX, DRZ, FW, VERBOSE;
	
	int kmax;
	int **recpos1, **recpos=NULL, nxrec=0, nyrec=0, nzrec=0;
	int itr=1, itr1=0, itr2=0, recflag=0, i, j, k, ifw, n;
	int nxrec1, nxrec2, nyrec1, nyrec2, nzrec1, nzrec2;
	float xrec, yrec, zrec;
	/*char rec_file_sub[STRING_SIZE]; */ /* variable not in use*/
	char bufferstring[10], buffer[STRING_SIZE];
	FILE *fpr;

	
	if (MYID==0) {
		switch (READREC) { /* read receiver positions from file */
			case 1:
				fprintf(fp,"\n Reading receiver positions from file: \n\t%s\n",REC_FILE);

				fpr=fopen(REC_FILE,"r");

				if (fpr==NULL) err(" Receiver file could not be opened !");

				*ntr=0;

				/* counts the number of receivers in the receiver file */
				while (fgets(buffer, STRING_SIZE, fpr)) {
					sscanf(buffer,"%s",bufferstring);

					/*testbuff=sscanf(buffer,"%s",bufferstring);
					fprintf(fp," bufferstring : _%s_with test=_%i_\n",bufferstring,testbuff); */
					/* checks if the line contains a '%'character which indicates a comment line,
					and if the reading of a string was successful, which is not the case for an empty line*/
					if ((strchr(buffer,'%')==0) && (sscanf(buffer,"%s",bufferstring)==1)) ++(*ntr);
				}

				rewind(fpr);

				recpos1=imatrix(1,3,1,*ntr);

				for (itr=1; itr<=*ntr; itr++) {
					fscanf(fpr,"%f%f%f\n",&xrec, &yrec, &zrec);
					/* "y" is used for the vertical coordinate*/
					recpos1[1][itr]=iround((xrec+REFREC[1])/DX);
					recpos1[2][itr]=iround((yrec+REFREC[2])/DY);
					recpos1[3][itr]=iround((zrec+REFREC[3])/DZ);
				}

				fclose(fpr);
				fprintf(fp," Message from function receiver (written by PE %d):\n",MYID);
				fprintf(fp," Number of receiver positions found: %i\n",*ntr);

				/* check if more than one receiver is located
								         at the same gridpoint */
				for (itr=1; itr<=(*ntr-1); itr++)
					for (itr1=itr+1; itr1<=*ntr; itr1++)
						if ((recpos1[1][itr]==recpos1[1][itr1])
						        && (recpos1[2][itr]==recpos1[2][itr1])
						        && (recpos1[3][itr]==recpos1[3][itr1]))
							recpos1[1][itr1]=-(++recflag);

				recpos=imatrix(1,3,1,*ntr-recflag);

				for (itr=1; itr<=*ntr; itr++)
					if (recpos1[1][itr]>0) {
						recpos[1][++itr2]=recpos1[1][itr];
						recpos[2][itr2]=recpos1[2][itr];
						recpos[3][itr2]=recpos1[3][itr];
					}

				*ntr=itr2;

				if ((recflag>0)||(itr2<(itr-1))) {
					fprintf(fp,"\n\n");
					fprintf(fp," Warning:\n");
					fprintf(fp," Several receivers located at the same gridpoint !\n");
					fprintf(fp," Number of receivers reduced to %i\n", *ntr);
					fprintf(fp,"\n\n");
				}

				free_imatrix(recpos1,1,3,1,*ntr);
				break;

			case 2 : /* REC ARRAY */
				ifw=FW+10;  /* frame width in gridpoints */

				if (BOUNDARY==1) ifw=0;

				*ntr=((1+(NZG-2*ifw)/DRZ)*(1+(NXG-2*ifw)/DRX))*REC_ARRAY;
				recpos=imatrix(1,3,1,*ntr);
				itr=0;

				for (n=0; n<=REC_ARRAY-1; n++) {
					j=iround((REC_ARRAY_DEPTH+REC_ARRAY_DIST*(float)n)/DY);

					for (k=ifw; k<=NZG-ifw; k+=DRZ)
						for (i=ifw; i<=NXG-ifw; i+=DRX) {
							itr++;
							recpos[1][itr]=i;
							recpos[2][itr]=j;
							recpos[3][itr]=k;
						}
				}

				break;

			case 0 :        /* straight horizontal or vertical line of receivers */
				
				if ((XREC1>XREC2) || ((YREC1>YREC2) ||(ZREC1>ZREC2))) {
					fprintf(fp," Coordinates of first receiver specified in input file :\n");
					fprintf(fp,"    %5.2f (x) , %5.2f (y) , %5.2f (z) :\n", XREC1,YREC1,ZREC1);
					fprintf(fp," Coordinates of last receiver specified in input file :\n");
					fprintf(fp,"    %5.2f (x) , %5.2f (y) , %5.2f (z) :\n", XREC2,YREC2,ZREC2);
					err("\n\n Receiver coordinates of first receiver should be equal/smaller than last receiver coordinates!");
				}

				nxrec1=iround(XREC1/DX); /* (nxrec1,nyrec1,nzrec1) and (nxrec2,nyrec2,nzrec2) */
				nyrec1=iround(YREC1/DY); /* are the positions of the first and last receiver*/
				nxrec2=iround(XREC2/DX); /* in gridpoints */
				nyrec2=iround(YREC2/DY);
				nzrec1=iround(ZREC1/DZ);
				nzrec2=iround(ZREC2/DZ);
				
				if ((abs(nyrec2-nyrec1)<=abs(nxrec2-nxrec1))||
				        (abs(nyrec2-nyrec1)<=abs(nzrec2-nzrec1))) {
					
					if (abs(nzrec2-nzrec1)<=abs(nxrec2-nxrec1)) {
						/* geophone-array horizontal x-dirextion */
						*ntr=iround((nxrec2-nxrec1)/NGEOPH)+1;
						recpos=imatrix(1,3,1,*ntr);
						for (nxrec=nxrec1; nxrec<=nxrec2; nxrec+=NGEOPH) {
							if ((nyrec2-nyrec1==0) && (nzrec2-nzrec1==0)) {
							nyrec=nyrec1;
							nzrec=nzrec1;	
							} else {	
							nyrec=nyrec1+((nyrec2-nyrec1)/(nxrec2-nxrec1)*(nxrec-nxrec1));
							nzrec=nzrec1+((nzrec2-nzrec1)/(nxrec2-nxrec1)*(nxrec-nxrec1));
							}
							itr=iround((nxrec-nxrec1)/NGEOPH)+1;
							
							recpos[1][itr]=nxrec;
							recpos[2][itr]=nyrec;
							recpos[3][itr]=nzrec;
						}
						
					} else { /* geophone-array horizontal z-direction */
						
						*ntr=iround((nzrec2-nzrec1)/NGEOPH)+1;
						/*fprintf(fp,"ntr=%d,nzrec2=%d, nzrec1=%d",ntr,nzrec2, nzrec1);*/
						recpos=imatrix(1,3,1,*ntr);

						for (nzrec=nzrec1; nzrec<=nzrec2; nzrec+=NGEOPH) {
							nyrec=nyrec1+((nyrec2-nyrec1)/(nzrec2-nzrec1)*(nzrec-nzrec1));
							nxrec=nxrec1+((nxrec2-nxrec1)/(nzrec2-nzrec1)*(nzrec-nzrec1));
							itr=iround((nzrec-nzrec1)/NGEOPH)+1;
							recpos[1][itr]=nxrec;
							recpos[2][itr]=nyrec;
							recpos[3][itr]=nzrec;
						}
					}

				} else {      /* receiver-line vertical */
					*ntr=iround((nyrec2-nyrec1)/NGEOPH)+1;
					recpos=imatrix(1,3,1,*ntr);
					
					for (nyrec=nyrec1; nyrec<=nyrec2; nyrec+=NGEOPH) {
						nxrec=nxrec1+((nxrec2-nxrec1)/(nyrec2-nyrec1)*(nyrec-nyrec1));
						nzrec=nzrec1+((nzrec2-nzrec1)/(nyrec2-nyrec1)*(nyrec-nyrec1));
						itr=iround((nyrec-nyrec1)/NGEOPH)+1;
						recpos[1][itr]=nxrec;
						recpos[2][itr]=nyrec;
						recpos[3][itr]=nzrec;
					}
				}


				break;
			default :
				err(" invalid READREC in receiver.c!");
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(ntr,1,MPI_INT,0,MPI_COMM_WORLD);

	if (MYID!=0) recpos=imatrix(1,3,1,*ntr);

	MPI_Bcast(&recpos[1][1],(*ntr)*3,MPI_INT,0,MPI_COMM_WORLD);

	if (MYID==0) {
		fprintf(fp,"\n **Message from function receiver (written by PE %d):\n",MYID);
		fprintf(fp," Number of receiver positions found: %i\n",*ntr);
		fprintf(fp," Receiver positions (in gridpoints) in the global model-system:\n");
		fprintf(fp," x  \ty \tz \n");
		fprintf(fp," -  \t- \t- \n");
		
		if (!VERBOSE) kmax=30; else kmax=*ntr;
		for (k=1; k<=kmax; k++)
			fprintf(fp," %5.2f   %5.2f   %5.2f\n",recpos[1][k]*DX,recpos[2][k]*DY,recpos[3][k]*DZ);
		if (!VERBOSE) fprintf(fp," ...\n Only %d of %d receiver positions displayed here",kmax,*ntr);
		fprintf(fp,"\n\n");
	}


	return recpos;
}
