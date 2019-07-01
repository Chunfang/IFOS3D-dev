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
 *   Write seismograms to disk
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "segy.h"

/* ****************************  UNDER CONSTRUCTION !!!  ***************************** */


void  outseis(FILE *fp, FILE *fpdata, int comp, float **section,
              int **recpos, int **recpos_loc, int ntr, float **srcpos,
              int nsrc, int ns, int seis_form) {
	/* declaration of extern variables */
	extern int NDT;
	extern float  DX, DY, DZ, TIME, DT, REFREC[4];

	extern int LITTLEBIG;

	/* declaration of local variables */
	int i,j, * pint;
	segy tr;
	int tracl ;
	float xr, yr, zr, y, z, scalefac, tfloat;  /*x*/
	float XS=0.0, YS=0.0, ZS=0.0;
	const float scale=3.0;
	

	if (nsrc==1) {
		/* only if one source position is specified in SOURCE_FILE,
			source coordinates are written into trace header fields */
		XS=srcpos[1][1];
		YS=srcpos[2][1];
		ZS=srcpos[3][1];
	}

	scalefac=pow(10.0,scale);

	switch (seis_form) {
	case 0 :
	case 1 : /* SU ~ (IEEE)  */

		/* fprintf(stderr,"BEEN HERE !!!\n"); */
		for (tracl=1; tracl<=ntr; tracl++) {
			xr=recpos[1][recpos_loc[4][tracl]]*DX;
			yr=recpos[2][recpos_loc[4][tracl]]*DY;
			zr=recpos[3][recpos_loc[4][tracl]]*DZ;
			/*x=xr-REFREC[1];*/
			y=yr-REFREC[2];
			z=zr-REFREC[3];
			tr.tracl=recpos_loc[4][tracl];      /* trace sequence number within line */
			tr.ep=comp;
			//tr.cdp=recpos_loc[4][tracl];
			tr.trid=1;           /* trace identification code: 1=seismic*/
			tr.offset=iround(sqrt((XS-xr)*(XS-xr)
			                      +(YS-yr)*(YS-yr)
			                      +(ZS-zr)*(ZS-zr))*scalefac);
			tr.gelev=iround(yr*scalefac);
			tr.sdepth=iround(YS*scalefac);   /* source depth (positive) */

			/* angle between receiver position and reference point
			   (sperical coordinate system: used for tunnel geometry) */
			tr.gdel=iround(atan2(-y,z)*180.0*scalefac/PI);
			tr.gwdep=iround(sqrt(z*z+y*y)*scalefac);

			tr.scalel=(short)(-scale);
			tr.scalco=(short)(-scale);
			tr.sx=iround(XS*scalefac);  /* X source coordinate */
			tr.sy=iround(ZS*scalefac);  /* Z source coordinate */

			/* group coordinates */
			tr.gx=iround(xr*scalefac);
			tr.gy=iround(zr*scalefac);


			tr.ns=(unsigned short)ns; /* number of samples in this trace */
			tr.dt=(unsigned short)iround(((float)NDT*DT)*1.0e6); /* sample interval in micro-seconds */
			tr.d1=(float)(TIME/ns);        /* sample spacing for non-seismic data */
			tr.d1=0;
			tr.tracr=0	;	/* trace sequence number within reel */

			tr.fldr=0       ;	/* field record number */

			tr.tracf=0      ;	/* trace number within field record */

			tr.ep=0         ;	/* energy source point number */

			tr.cdpt=0       ;	/* trace number within CDP ensemble */


			tr.nvs=0       	;   /* number of vertically summed traces (see vscode
			   			in bhed structure) */

			tr.nhs=0       	;   /* number of horizontally summed traces (see vscode
			  			in bhed structure) */

			tr.duse=0     	;   /* data use:
						1 = production
						2 = test */

			tr.gdel=0      	; /* datum elevation at receiver group */

			tr.sdel=0      	; /* datum elevation at source */

			tr.gwdep=0     	; /* water depth at receiver group */

			tr.counit=0    	;   /* coordinate units code:
						for previous four entries
						1 = length (meters or feet)
						2 = seconds of arc (in this case, the
						X values are longitude and the Y values
						are latitude, a positive value designates
						the number of seconds east of Greenwich
						or north of the equator */

			tr.wevel=0     ;	/* weathering velocity */

			tr.swevel=0    ;	/* subweathering velocity */

			tr.sut=0       ;	/* uphole time at source */

			tr.gut=0       ;	/* uphole time at receiver group */

			tr.sstat=0     ;	/* source static correction */

			tr.gstat=0     ;	/* group static correction */

			tr.tstat=0     ;	/* total static applied */

			tr.laga=0      ; /* lag time A, time in ms between end of 240-
			   			byte trace identification header and time
			   			break, positive if time break occurs after
			   			end of header, time break is defined as
			   			the initiation pulse which maybe recorded
			   			on an auxiliary trace or as otherwise
			   			specified by the recording system */

			tr.lagb=0     	; /* lag time B, time in ms between the time break
			   			and the initiation time of the energy source,
			   			may be positive or negative */

			tr.delrt=0     	; /* delay recording time, time in ms between
			   			initiation time of energy source and time
			   			when recording of data samples begins
			   			(for deep water work if recording does not
			   			start at zero time) */

			tr.muts=0      ; /* mute time--start */

			tr.mute=0      ; /* mute time--end */

			tr.gain=0      ; /* gain type of field instruments code:
						1 = fixed
						2 = binary
						3 = floating point
						4 ---- N = optional use */

			tr.igc=0       ; /* instrument gain constant */

			tr.igi=0       ; /* instrument early or initial gain */

			tr.corr=0      ; /* correlated:
						1 = no
						2 = yes */

			tr.sfs=0       ; /* sweep frequency at start */

			tr.sfe=0       ; /* sweep frequency at end */

			tr.slen=0      ; /* sweep length in ms */

			tr.styp=0      ; /* sweep type code:
						1 = linear
						2 = cos-squared
						3 = other */

			tr.stas=0      	; /* sweep trace length at start in ms */

			tr.stae=0      	; /* sweep trace length at end in ms */

			tr.tatyp=0     	; /* taper type: 1=linear, 2=cos^2, 3=other */

			tr.afilf=0     	; /* alias filter frequency if used */

			tr.afils=0     	; /* alias filter slope */

			tr.nofilf=0    	; /* notch filter frequency if used */

			tr.nofils=0   	; /* notch filter slope */

			tr.lcf=0      	; /* low cut frequency if used */

			tr.hcf=0      	; /* high cut frequncy if used */

			tr.lcs=0       	; /* low cut slope */

			tr.hcs=0       	; /* high cut slope */

			tr.year=0      	; /* year data recorded */

			tr.day=0       	; /* day of year */

			tr.hour=0     	; /* hour of day (24 hour clock) */

			tr.minute=0    	; /* minute of hour */

			tr.sec=0       	; /* second of minute */

			tr.timbas=0    	; /* time basis code:
						1 = local
						2 = GMT
						3 = other */

			tr.trwf=0      	; /* trace weighting factor, defined as 1/2^N
			   			volts for the least sigificant bit */

			tr.grnors=0   	; /* geophone group number of roll switch
			   			position one */

			tr.grnofr=0    	; /* geophone group number of trace one within
			   			original field record */

			tr.grnlof=0    	; /* geophone group number of last trace within
			   			original field record */

			tr.gaps=0      	;  /* gap size (total number of groups dropped) */

			tr.otrav=0     	;  /* overtravel taper code:
						1 = down (or behind)
						2 = up (or ahead) */

			/* local assignments */

			tr.f1=0.0;	/* first sample location for non-seismic data */

			tr.d2=0.0;	/* sample spacing between traces */

			tr.f2=0.0;	/* first trace location */

			tr.ungpow=0.0;	/* negative of power used for dynamic
			   			range compression */

			tr.unscale=0.0;	/* reciprocal of scaling factor to normalize
			   			range */
			tr.ntr=0      ;   /* number of traces */

			tr.mark=0     ;

			for (j=1; j<=ns; j++) {
				tr.data[j]=section[tracl][j];
			}

			fwrite(&tr,240,1,fpdata);
			fwrite(&tr.data[1],4,ns,fpdata);
		}

		break;

	case 2 :


		for (j=1; j<=ns; j++) {     /*ASCII ONE COLUMN PER TRACE */
			for (i=1; i<=ntr; i++) {
				fprintf(fpdata,"%e\t", section[i][j]);
			}

			fprintf(fpdata,"\n");
		}



		break;

	case 3 :                             /*BINARY */
		if (!LITTLEBIG) { /* OUTPUT NATIVE FLOATS */
			for (i=1; i<=ntr; i++) for (j=1; j<=ns; j++) {
					fwrite(&section[i][j],sizeof(float),1,fpdata);
				}

		} else {
				/* SWAP FLOATS */
				for (i=1; i<=ntr; i++) for (j=1; j<=ns; j++) {
						tfloat=section[i][j];
						pint=(int *) &tfloat;
						*pint=((*pint>>24)&0xff)|((*pint&0xff)<<24)|((*pint>>8)&0xff00)|((*pint&0xff00)<<8);
							/* explanation:
							* (*pint>>24)&0xff shifts 3 bytes (3*8) to the right and takes with &=and only the least significant byte (other bytes 0) 0 0 0 1
							(*pint&0xff)<<24 sets the first 3 bytes to zero (takes only the least significant byte) than shifts 3 bytes to the left (0 0 0 1 -> 1 0 0 0) 
							(*pint>>8)&0xff00 shifts 1 byte to the right and than takes only 3rd byte (0 0 1 0)
							(*pint&0xff00)<<8 sets bytes 1,2,4 to zero and shifts 1 byte to the left (0 1 0 0)
							The |=or operator adds all for parts to the swaped float
							*/
						fwrite(&tfloat,sizeof(float),1,fpdata);
					}

		}

		break;

	default :
		fprintf(fp," Don't know data format for seismograms ! Choose SEIS_FORMAT 1-3 \n");
		fprintf(fp," No output written. ");
	}

	fclose(fpdata);


}
