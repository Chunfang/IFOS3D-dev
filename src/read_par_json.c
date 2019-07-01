/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 *
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 *
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
 -----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   program IFOS, reading input-parameters from input-file or stdin
 *  ----------------------------------------------------------------------*/

#include <unistd.h>
#include "fd.h"

char **varname_list,* *value_list;

void read_par_json(FILE *fp, char *fileinp) {

	/* declaration of extern variables */
	extern int   NX, NY, NZ, SOURCE_SHAPE, SOURCE_TYPE, SNAP, SNAP_FORMAT, SNAP_PLANE;
	extern int DRX, DRZ, L, SRCREC, FDORDER, FW, FDCOEFF;
	extern float DX, DY, DZ, TIME, DT, TS, *FL, TAU, PLANE_WAVE_DEPTH, PHI;
	extern float XREC1, XREC2, YREC1, YREC2, ZREC1, ZREC2, ALPHA, BETA;
	extern float REC_ARRAY_DEPTH, REC_ARRAY_DIST;
	extern int SEISMO, NDT, NDTSHIFT, NGEOPH, SEIS_FORMAT, FREE_SURF, READMOD, MOD_FORMAT, READREC, RUN_MULTIPLE_SHOTS;
	extern int BOUNDARY, REC_ARRAY, IDX, IDY, IDZ, ABS_TYPE;
	extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4], DAMPING, FPML, VPPML;
	extern char  MFILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE];
	extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char SEIS_FILE[STRING_SIZE],GRAD_FILE[STRING_SIZE], SEIS_OBS_FILE[STRING_SIZE],INV_FILE[STRING_SIZE];
	extern int NPROCX,NPROCY,NPROCZ;
	extern int LITTLEBIG;

	extern float REFSRC[3], SRCTSHIFT;
	extern int SRC_MF, SIGNAL_FORMAT;
	extern char MOD_OUT_FILE[STRING_SIZE], HESS_FILE[STRING_SIZE];
	extern int METHOD;
	extern int ITMIN, ITMAX, FILT, NFMAX, TAST, NSHOTS_STEP, DAMPTYPE, HESS, READ_HESS, REC_HESS,EXTOBS,LBFGS;
	/*extern float F_INV;*/
	extern float TESTSTEP, WATER_HESS[3], WEIGHT[3], VP0, VS0, RHO0;
	extern int BFGSNUM, NUMPAR;
	extern int MYID;
	extern int VERBOSE;
	/* definition of local variables */

	int number_readobjects=0,fserr=0;
	char errormessage[STRING_SIZE2];

	char **varname_list, ** value_list;

	if (MYID == 0) {

		/* allocate first object in list */
		varname_list = malloc(STRING_SIZE2*sizeof(char *));
		value_list = malloc(STRING_SIZE2*sizeof(char *));

		/* read in objects from file */
		number_readobjects=read_objects_from_intputfile(fp, fileinp, varname_list, value_list);
		fprintf(fp,"\n From input file %s, %i objects have been read in. \n",fileinp, number_readobjects);

		/* print objects to screen */
		fprintf(fp, "\n ===========================================================");
		fprintf(fp, "\n =   List of Parameters read by the built in Json Parser   =");
		print_objectlist_screen(fp, number_readobjects, varname_list, value_list);

		/* extract variables form object list */

		/*=================================
		 section general grid and discretization parameters
		 =================================*/
		if (get_int_from_objectlist("NPROCX",number_readobjects,&NPROCX,varname_list, value_list)) {
			err("Variable NPROCX could not be retrieved from the json input file!");
		}

		if (get_int_from_objectlist("NPROCY",number_readobjects,&NPROCY,varname_list, value_list)) {
			err("Variable NPROCY could not be retrieved from the json input file!");
		}

		if (get_int_from_objectlist("NPROCZ",number_readobjects,&NPROCZ,varname_list, value_list)) {
			err("Variable NPROCZ could not be retrieved from the json input file!");
		}


		if (get_int_from_objectlist("NX",number_readobjects,&NX,varname_list, value_list)) {
			err("Variable NX could not be retrieved from the json input file!");
		}

		if (get_int_from_objectlist("NY",number_readobjects,&NY,varname_list, value_list)) {
			err("Variable NY could not be retrieved from the json input file!");
		}

		if (get_int_from_objectlist("NZ",number_readobjects,&NZ,varname_list, value_list)) {
			err("Variable NZ could not be retrieved from the json input file!");
		}

		if (get_float_from_objectlist("DX",number_readobjects,&DX,varname_list, value_list)) {
			err("Variable DX could not be retrieved from the json input file!");
		}

		if (get_float_from_objectlist("DY",number_readobjects,&DY,varname_list, value_list)) {
			err("Variable DY could not be retrieved from the json input file!");
		}

		if (get_float_from_objectlist("DZ",number_readobjects,&DZ,varname_list, value_list)) {
			err("Variable DZ could not be retrieved from the json input file!");
		}


		if (get_int_from_objectlist("FDORDER",number_readobjects,&FDORDER,varname_list, value_list)) {
			err("Variable FDORDER could not be retrieved from the json input file!");
		}

		if (get_int_from_objectlist("FDCOEFF",number_readobjects,&FDCOEFF,varname_list, value_list)) {
			err("Variable FDCOEFF could not be retrieved from the json input file!");
		}


		if (get_float_from_objectlist("TIME",number_readobjects,&TIME,varname_list, value_list)) {
			err("Variable TIME could not be retrieved from the json input file!");
		}

		if (get_float_from_objectlist("DT",number_readobjects,&DT,varname_list, value_list)) {
			err("Variable DT could not be retrieved from the json input file!");
		}

		/*=================================
		 section source parameters
		 =================================*/
		fprintf(fp," The following default values are set:\n");
		fprintf(fp," =====================================\n\n");


		if (get_int_from_objectlist("SOURCE_SHAPE",number_readobjects,&SOURCE_SHAPE,varname_list, value_list)) {
			err(" Variable SOURCE_SHAPE could not be retrieved from the json input file!");
		}

		else {
			if (SOURCE_SHAPE==3) {
				if (get_string_from_objectlist("SIGNAL_FILE",number_readobjects,SIGNAL_FILE,varname_list, value_list)) {
					err(" Variable SIGNAL_FILE could not be retrieved from the json input file!");

				} else {
					if (get_int_from_objectlist("SIGNAL_FORMAT",number_readobjects,&SIGNAL_FORMAT,varname_list, value_list)) {
						SIGNAL_FORMAT=1;
					}
				}
			}
		}

		if (get_int_from_objectlist("SOURCE_TYPE",number_readobjects,&SOURCE_TYPE,varname_list, value_list)) {
			err(" Variable SOURCE_TYPE could not be retrieved from the json input file!");
		}

		else {
			if (SOURCE_TYPE==5) {
				if (get_float_from_objectlist("ALPHA",number_readobjects,&ALPHA,varname_list, value_list)) {
					err(" Variable ALPHA could not be retrieved from the json input file!");
				}

				if (get_float_from_objectlist("BETA",number_readobjects,&BETA,varname_list, value_list)) {
					err(" Variable BETA could not be retrieved from the json input file!");
				}
			}

		}

		if (get_int_from_objectlist("SRCREC",number_readobjects,&SRCREC,varname_list, value_list)) {
			err(" Variable SRCREC could not be retrieved from the json input file!");

		} else {
			if (get_int_from_objectlist("SRC_MF",number_readobjects,&SRC_MF,varname_list, value_list)) {
				SRC_MF=0;
			}

			if (get_float_from_objectlist("REFSRCX",number_readobjects,&REFSRC[0],varname_list, value_list)) {
				REFSRC[0]=0.0;
				fprintf(fp," Variable REFSRCX is set to default value %.1f.\n",REFSRC[0]);
			}

			if (get_float_from_objectlist("REFSRCY",number_readobjects,&REFSRC[1],varname_list, value_list)) {
				REFSRC[1]=0.0;
				fprintf(fp," Variable REFSRCY is set to default value %.1f.\n",REFSRC[1]);
			}

			if (get_float_from_objectlist("REFSRCZ",number_readobjects,&REFSRC[2],varname_list, value_list)) {
				REFSRC[2]=0.0;
				fprintf(fp," Variable REFSRCZ is set to default value %.1f.\n",REFSRC[2]);
			}


			if (SRCREC==1) {
				if (get_string_from_objectlist("SOURCE_FILE",number_readobjects,SOURCE_FILE,varname_list, value_list)) {
					err(" Variable SOURCE_FILE could not be retrieved from the json input file!");
				}

				if (get_int_from_objectlist("RUN_MULTIPLE_SHOTS",number_readobjects,&RUN_MULTIPLE_SHOTS,varname_list, value_list)) {
					err(" Variable RUN_MULTIPLE_SHOTS could not be retrieved from the json input file!");

				} else {
					if (get_float_from_objectlist("SRCTSHIFT",number_readobjects,&SRCTSHIFT,varname_list, value_list)) {
						SRCTSHIFT=0.0;
					}

				}
			}


			if (SRCREC==2) {
				if (get_float_from_objectlist("PLANE_WAVE_DEPTH",number_readobjects,&PLANE_WAVE_DEPTH,varname_list, value_list)) {
					err(" Variable PLANE_WAVE_DEPTH could not be retrieved from the json input file!");

				} else {
					if (PLANE_WAVE_DEPTH>0) {
						if (get_float_from_objectlist("PHI",number_readobjects,&PHI,varname_list, value_list)) {
							err("Variable PHI could not be retrieved from the json input file!");
						}

						if (get_float_from_objectlist("TS",number_readobjects,&TS,varname_list, value_list)) {
							err("Variable TS could not be retrieved from the json input file!");

						}
					}
				}
			}
		} /* end of SRCREC */



		/*=================================
		 section general model and log parameters
		 =================================*/
		if (get_int_from_objectlist("VERBOSE",number_readobjects,&VERBOSE,varname_list, value_list)) {
			VERBOSE=1;
			fprintf(fp,"Variable VERBOSE is set to value %d.\n",VERBOSE);
		}

		if (get_int_from_objectlist("READMOD",number_readobjects,&READMOD,varname_list, value_list)) {
			err("Variable READMOD could not be retrieved from the json input file!");

		} else {
			if (get_int_from_objectlist("MOD_FORMAT",number_readobjects,&MOD_FORMAT,varname_list, value_list)) {
				MOD_FORMAT=0;
			}

			if (get_string_from_objectlist("MFILE",number_readobjects,MFILE,varname_list, value_list)) {
				err("Variable MFILE could not be retrieved from the json input file!");
			}
		}


		if (get_int_from_objectlist("L",number_readobjects,&L,varname_list, value_list)) {
			L=0;
			fprintf(fp," Variable L is set to default value %d.\n",L);

		} else {
			FL=vector(1,L);

			switch (L) {
				case 0:
					break;

				case 1:
					if (get_float_from_objectlist("FL1",number_readobjects,&FL[1],varname_list, value_list)) {
						err("Variable FL1 could not be retrieved from the json input file!");
					}

					break;

				default:
					err("More than four relaxation Parameter (L>1) are not implemented yet!");
					break;
			}

			if (L) {
				if (get_float_from_objectlist("TAU",number_readobjects,&TAU,varname_list, value_list)) {
					err("Variable TAU could not be retrieved from the json input file!");
				}
			}
		}

		/*=================================
			  section boundary parameters
		  =================================*/

		if (get_int_from_objectlist("FREE_SURF",number_readobjects,&FREE_SURF,varname_list, value_list)) {
			err("Variable FREE_SURF could not be retrieved from the json input file!");
		}

		if (get_int_from_objectlist("BOUNDARY",number_readobjects,&BOUNDARY,varname_list, value_list)) {
			BOUNDARY=0;
			fprintf(fp," Variable BOUNDARY is set to default value %d.\n",BOUNDARY);
		}

		if (get_int_from_objectlist("ABS_TYPE",number_readobjects,&ABS_TYPE,varname_list, value_list)) {
			err("Variable ABS_TYPE could not be retrieved from the json input file!");
		}

		if (ABS_TYPE==1) {
			if (get_float_from_objectlist("FPML",number_readobjects,&FPML,varname_list, value_list)) {
				err("Variable FPML could not be retrieved from the json input file!");
			}

			if (get_float_from_objectlist("VPPML",number_readobjects,&VPPML,varname_list, value_list)) {
				err("Variable VPPML could not be retrieved from the json input file!");
			}
		}

		if (get_int_from_objectlist("FW",number_readobjects,&FW,varname_list, value_list)) {
			err("Variable FW could not be retrieved from the json input file!");
		}

		if (ABS_TYPE==2) {
			if (get_float_from_objectlist("DAMPING",number_readobjects,&DAMPING,varname_list, value_list)) {
				err("Variable DAMPING could not be retrieved from the json input file!");
			}
		}


		/*=================================
		section snapshot parameters
		=================================*/
		if (get_int_from_objectlist("SNAP",number_readobjects,&SNAP,varname_list, value_list)) {
			err("Variable SNAP not be retrieved from the json input file!");

		} else {
			if (SNAP>0) {
				if (get_int_from_objectlist("SNAP_FORMAT",number_readobjects,&SNAP_FORMAT,varname_list, value_list)) {
					err("Variable SNAP_FORMAT could not be retrieved from the json input file!");
				}

				if (get_float_from_objectlist("TSNAP1",number_readobjects,&TSNAP1,varname_list, value_list)) {
					err("Variable TSNAP1 could not be retrieved from the json input file!");
				}

				if (get_float_from_objectlist("TSNAP2",number_readobjects,&TSNAP2,varname_list, value_list)) {
					err("Variable TSNAP2 could not be retrieved from the json input file!");
				}

				if (get_float_from_objectlist("TSNAPINC",number_readobjects,&TSNAPINC,varname_list, value_list)) {
					err("Variable TSNAPINC could not be retrieved from the json input file!");
				}

				if (get_string_from_objectlist("SNAP_FILE",number_readobjects,SNAP_FILE,varname_list, value_list)) {
					err("Variable SNAP_FILE could not be retrieved from the json input file!");
				}

			}
		}

		if (SNAP==3) {
			if (get_int_from_objectlist("SNAP_PLANE",number_readobjects,&SNAP_PLANE,varname_list, value_list)) {
				err("Variable SNAP_PLANE could not be retrieved from the json input file!");
			}
		}

		/* increments are read in any case, because they will be also used as increment for model output */
		if (get_int_from_objectlist("IDX",number_readobjects,&IDX,varname_list, value_list)) {
			IDX=1;
			fprintf(fp," Variable IDX is set to default value %d.\n",IDX);
		}

		if (get_int_from_objectlist("IDY",number_readobjects,&IDY,varname_list, value_list)) {
			IDY=1;
			fprintf(fp," Variable IDY is set to default value %d.\n",IDY);
		}

		if (get_int_from_objectlist("IDZ",number_readobjects,&IDZ,varname_list, value_list)) {
			IDZ=1;
			fprintf(fp," Variable IDZ is set to default value %d.\n",IDZ);
		}

		/*=================================
		 section seismogramm parameters
		 =================================*/
		if (get_int_from_objectlist("SEISMO",number_readobjects,&SEISMO,varname_list, value_list)) {
			err("Variable SEISMO could not be retrieved from the json input file!");
		}

		else {
			if (SEISMO>0) {
				if (get_string_from_objectlist("SEIS_FILE",number_readobjects,SEIS_FILE,varname_list, value_list)) {
					err("Variable SEIS_FILE could not be retrieved from the json input file!");
				}

				if (get_int_from_objectlist("READREC",number_readobjects,&READREC,varname_list, value_list)) {
					err("Variable READREC could not be retrieved from the json input file!");
				}

				else {
					switch (READREC) {
						case 0 : /*Receiver line*/
							if (get_float_from_objectlist("XREC1",number_readobjects,&XREC1,varname_list, value_list)) {
								err("Variable XREC1 could not be retrieved from the json input file!");
							}

							if (get_float_from_objectlist("XREC2",number_readobjects,&XREC2,varname_list, value_list)) {
								err("Variable XREC2T could not be retrieved from the json input file!");
							}

							if (get_float_from_objectlist("YREC1",number_readobjects,&YREC1,varname_list, value_list)) {
								err("Variable YREC1 could not be retrieved from the json input file!");
							}

							if (get_float_from_objectlist("YREC2",number_readobjects,&YREC2,varname_list, value_list)) {
								err("Variable YREC2 could not be retrieved from the json input file!");
							}

							if (get_float_from_objectlist("ZREC1",number_readobjects,&ZREC1,varname_list, value_list)) {
								err("Variable ZREC1 could not be retrieved from the json input file!");
							}

							if (get_float_from_objectlist("ZREC2",number_readobjects,&ZREC2,varname_list, value_list)) {
								err("Variable ZREC2 could not be retrieved from the json input file!");
							}


							if (get_int_from_objectlist("NGEOPH",number_readobjects,&NGEOPH,varname_list, value_list)) {
								err("Variable NGEOPH could not be retrieved from the json input file!");
							}

							break;

						case 1 : /*Receiver from file*/
							if (get_string_from_objectlist("REC_FILE",number_readobjects,REC_FILE,varname_list, value_list)) {
								err("Variable REC_FILE could not be retrieved from the json input file!");
							}

							if (get_float_from_objectlist("REFRECX",number_readobjects,&REFREC[1],varname_list, value_list)) {
								REFREC[1]=0.0;
								fprintf(fp," Variable REFREC is set to default value %.1f.\n",REFREC[1]);
							}

							if (get_float_from_objectlist("REFRECY",number_readobjects,&REFREC[2],varname_list, value_list)) {
								REFREC[2]=0.0;
								fprintf(fp," Variable REFREC is set to default value %.1f.\n",REFREC[1]);
							}

							if (get_float_from_objectlist("REFRECZ",number_readobjects,&REFREC[3],varname_list, value_list)) {
								REFREC[3]=0.0;
								fprintf(fp," Variable REFREC is set to default value %.1f.\n",REFREC[1]);
							}

							break;

						case 2: /*Receiver array*/
							if (get_int_from_objectlist("REC_ARRAY",number_readobjects,&REC_ARRAY,varname_list, value_list)) {
								err("Variable REC_ARRAY could not be retrieved from the json input file!");
							}

							if (get_float_from_objectlist("REC_ARRAY_DEPTH",number_readobjects,&REC_ARRAY_DEPTH,varname_list, value_list)) {
								err("Variable REC_ARRAY_DEPTH could not be retrieved from the json input file!");
							}

							if (get_float_from_objectlist("REC_ARRAY_DIST",number_readobjects,&REC_ARRAY_DIST,varname_list, value_list)) {
								err("Variable REC_ARRAY_DIST could not be retrieved from the json input file!");
							}

							if (get_int_from_objectlist("DRX",number_readobjects,&DRX,varname_list, value_list)) {
								err("Variable DRX could not be retrieved from the json input file!");
							}

							if (get_int_from_objectlist("DRZ",number_readobjects,&DRZ,varname_list, value_list)) {
								err("Variable DRZ could not be retrieved from the json input file!");
							}

							break;

						default :
							err("Please choose READREC=0 (Receiver Line) ,READREC=1 (Receiver from file) or READREC=2(Receiver Array)");


					}
				}

				if (READREC!=1) {
					REFREC[1]=0.0;
					REFREC[2]=0.0;
					REFREC[3]=0.0;
					fprintf(fp," Variable REFREC is set to default value (%.1f,%.1f,%.1f).\n",REFREC[1],REFREC[2],REFREC[3]);
				}
			}


			/* --------output ----------
			 *------------------------*/

			if (get_int_from_objectlist("NDT",number_readobjects,&NDT,varname_list, value_list)) {
				NDT=1;
				fprintf(fp," Variable NDT is set to default value %d.\n",NDT);
			}

			if (get_int_from_objectlist("NDTSHIFT",number_readobjects,&NDTSHIFT,varname_list, value_list)) {
				NDTSHIFT=0;
				fprintf(fp," Variable NDTSHIFT is set to default value %d.\n",NDT);
			}

			if (get_int_from_objectlist("SEIS_FORMAT",number_readobjects,&SEIS_FORMAT,varname_list, value_list)) {
				err("Variable SEIS_FORMAT could not be retrieved from the json input file!");

			} 
		}/*end of seismo*/


		if (get_int_from_objectlist("LITTLEBIG",number_readobjects,&LITTLEBIG,varname_list, value_list)) {
			LITTLEBIG=0;
		}


		/*=================================
		 section inversion parameters
		 =================================*/


		if (get_int_from_objectlist("METHOD",number_readobjects,&METHOD,varname_list, value_list)) {
			err("Variable METHOD could not be retrieved from the json input file!");
		}

		else {
			if (METHOD==1) {	/* FWI is calculated */


				/*=================================
				section General
				=================================*/

				if (get_int_from_objectlist("ITMIN",number_readobjects,&ITMIN,varname_list, value_list)) {
					err("Variable ITMIN could not be retrieved from the json input file!");
				}

				if (get_int_from_objectlist("ITMAX",number_readobjects,&ITMAX,varname_list, value_list)) {
					err("Variable ITMAX could not be retrieved from the json input file!");
				}

				if (get_int_from_objectlist("FILT",number_readobjects,&FILT,varname_list, value_list)) {
					err("Variable FILT could not be retrieved from the json input file!");
				}


				if (get_int_from_objectlist("NFMAX",number_readobjects,&NFMAX,varname_list, value_list)) {
					err("Variable NFMAX could not be retrieved from the json input file!");
				}

				if (get_int_from_objectlist("TAST",number_readobjects,&TAST,varname_list, value_list)) {
					err("Variable TAST could not be retrieved from the json input file!");
				}

				if (get_float_from_objectlist("VP0",number_readobjects,&VP0,varname_list, value_list)) {
					err("Variable VP0 could not be retrieved from the json input file!");
				}

				if (get_float_from_objectlist("VS0",number_readobjects,&VS0,varname_list, value_list)) {
					err("Variable VS0 could not be retrieved from the json input file!");
				}

				if (get_float_from_objectlist("RHO0",number_readobjects,&RHO0,varname_list, value_list)) {
					err("Variable RHO0 could not be retrieved from the json input file!");
				}

				if (get_float_from_objectlist("RHO0",number_readobjects,&RHO0,varname_list, value_list)) {
					err("Variable RHO0 could not be retrieved from the json input file!");
				}

				if (get_float_from_objectlist("WEIGHT_VP",number_readobjects,&WEIGHT[0],varname_list, value_list)) {
					err("Variable WEIGHT_VP could not be retrieved from the json input file!");
				}

				if (get_float_from_objectlist("WEIGHT_VS",number_readobjects,&WEIGHT[1],varname_list, value_list)) {
					err("Variable WEIGHT_VS could not be retrieved from the json input file!");
				}

				if (get_float_from_objectlist("WEIGHT_RHO",number_readobjects,&WEIGHT[2],varname_list, value_list)) {
					err("Variable WEIGHT_RHO could not be retrieved from the json input file!");
				}



				/*=================================
				section Steplength estimation + Gradient preconditioning
				=================================*/

				if (get_int_from_objectlist("NSHOTS_STEP",number_readobjects,&NSHOTS_STEP,varname_list, value_list)) {
					err("Variable NSHOTS_STEP could not be retrieved from the json input file!");
				}

				if (get_float_from_objectlist("TESTSTEP",number_readobjects,&TESTSTEP,varname_list, value_list)) {
					err("Variable TESTSTEP could not be retrieved from the json input file!");
				}

				if (get_int_from_objectlist("DAMPTYPE",number_readobjects,&DAMPTYPE,varname_list, value_list)) {
					err("Variable DAMPTYPE could not be retrieved from the json input file!");
				}

				/*=================================
				section Hessian + L-BFGS
				=================================*/

				if (get_int_from_objectlist("HESS",number_readobjects,&HESS,varname_list, value_list)) {
					err("Variable HESS could not be retrieved from the json input file!");

				} else {
					if (HESS) {
						if (get_int_from_objectlist("READ_HESS",number_readobjects,&READ_HESS,varname_list, value_list)) {
							err("Variable READ_HESS could not be retrieved from the json input file!");
						}

						if (get_int_from_objectlist("REC_HESS",number_readobjects,&REC_HESS,varname_list, value_list)) {
							err("Variable REC_HESS could not be retrieved from the json input file!");
						}

						if (get_float_from_objectlist("WATER_HESS_VP",number_readobjects,&WATER_HESS[0],varname_list, value_list)) {
							err("Variable WATER_HESS_VP could not be retrieved from the json input file!");
						}

						if (get_float_from_objectlist("WATER_HESS_VS",number_readobjects,&WATER_HESS[1],varname_list, value_list)) {
							err("Variable WATER_HESS_VS could not be retrieved from the json input file!");
						}

						if (get_float_from_objectlist("WATER_HESS_RHO",number_readobjects,&WATER_HESS[2],varname_list, value_list)) {
							err("Variable WATER_HESS_RHO could not be retrieved from the json input file!");
						}
					}
				}

				if (get_int_from_objectlist("LBFGS",number_readobjects,&LBFGS,varname_list, value_list)) {
					err("Variable LBFGS could not be retrieved from the json input file!");

				} else {
					if (LBFGS) {

						if (get_int_from_objectlist("NUMPAR",number_readobjects,&NUMPAR,varname_list, value_list)) {
							err("Variable NUMPAR could not be retrieved from the json input file!");
						}

						if (get_int_from_objectlist("BFGSNUM",number_readobjects,&BFGSNUM,varname_list, value_list)) {
							err("Variable BFGSNUM could not be retrieved from the json input file!");
						}
					}
				}

				/*=================================
				section In- and Output Files
				=================================*/

				if (get_string_from_objectlist("GRAD_FILE",number_readobjects,GRAD_FILE,varname_list, value_list)) {
					err("Variable GRAD_FILE could not be retrieved from the json input file!");
				}

				if (get_string_from_objectlist("MOD_OUT_FILE",number_readobjects,MOD_OUT_FILE,varname_list, value_list)) {
					err("Variable MOD_OUT_FILE could not be retrieved from the json input file!");
				}

				if (get_string_from_objectlist("SEIS_OBS_FILE",number_readobjects,SEIS_OBS_FILE,varname_list, value_list)) {
					err("Variable SEIS_OBS_FILE could not be retrieved from the json input file!");
				}

				if (get_int_from_objectlist("EXTOBS",number_readobjects,&EXTOBS,varname_list, value_list)) {
					err("Variable EXTOBS could not be retrieved from the json input file!");
				}

				if (get_string_from_objectlist("INV_FILE",number_readobjects,INV_FILE,varname_list, value_list)) {
					err("Variable INV_FILE could not be retrieved from the json input file!");
				}

				if (HESS) {
					if (get_string_from_objectlist("HESS_FILE",number_readobjects,HESS_FILE,varname_list, value_list)) {
						err("Variable HESS_FILE could not be retrieved from the json input file!");
					}
				}
			} /* end if (METHOD==1) */

			else {/* only forward modeling is applied */

				ITMIN=1;
				fprintf(fp," Variable ITMIN is set to default value %d.\n",ITMIN);
				ITMAX=1;
				fprintf(fp," Variable ITMAX is set to default value %d.\n",ITMAX);

				if (get_int_from_objectlist("FILT",number_readobjects,&FILT,varname_list, value_list)) {
					FILT=0;
					fprintf(fp," Variable FILT is set to default value %d.\n",FILT);

				if (get_string_from_objectlist("MOD_OUT_FILE",number_readobjects,MOD_OUT_FILE,varname_list, value_list)) {
					err("Variable MOD_OUT_FILE could not be retrieved from the json input file!");
				}

				}

			}

		}

		fprintf(fp,"\n End of setting default values\n");
		fprintf(fp," =====================================\n\n");


		/********************************************/
		/* Check files and directories if necessary */
		/********************************************/

		/* signal file */
		if (SOURCE_SHAPE == 3) {
			if (access(SIGNAL_FILE,0) != 0) {
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The signal file does not exist!\n");
				fprintf(fp, "        File name: <%s>", SIGNAL_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;

			} else if (access(SIGNAL_FILE,4) != 0) {
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The signal file does not have read access!\n");
				fprintf(fp, "        File name: <%s>", SIGNAL_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
		}

		/* source file */
		if (SRCREC==1) {
			if (access(SOURCE_FILE,0) != 0) {
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The source file does not exist!\n");
				fprintf(fp, "        File name: <%s>", SOURCE_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;

			} else if (access(SOURCE_FILE,4) != 0) {
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The source file does not have read access!\n");
				fprintf(fp, "        File name: <%s>", SOURCE_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
		}


		/* receiver file */
		if (READREC==1) {
			if (access(REC_FILE,0) != 0) {
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The receiver file does not exist!\n");
				fprintf(fp, "        File name: <%s>", REC_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;

			} else if (access(REC_FILE,4) != 0) {
				fprintf(fp, "\n==================================================================\n");
				fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
				fprintf(fp, "        The receiver file does not have read access!\n");
				fprintf(fp, "        File name: <%s>", REC_FILE);
				fprintf(fp, "\n==================================================================\n");
				fserr = 1;
			}
		}


		/********************************************/
		/* ERROR                                    */
		/********************************************/
		if (fserr) {
			fprintf(fp, "\n");
			sprintf(errormessage, "\n  in: <read_par_json.c> \n");
			err(errormessage);
		}


	} /* End of if(MYID==0) */
}


