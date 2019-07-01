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
 *   loop over snapshotfiles which have to be merged.
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "globvar.h"      /* definition of global variables  */


int main(int argc, char **argv) {

	int nsnap;
	char *fileinp="";
	//FILE *FP;
	fileinp = argv[1];


	if ((FP=fopen(fileinp,"r"))==NULL) {
		err(" Opening input file failed.");

	} else {
		printf(" Opening input file was successful.\n\n");
	}

	/* read parameters from parameter-file */

	if (strstr(fileinp,".json")) {
		//read json formated input file
		read_par_json(stdout, fileinp);
		fclose(FP);

	} else {
		//read "old" input file *.inp, might not work in future
		err(" Old Input files (.inp) are no longer supported. \n Please use .json input files instead. \n\n");

	}



	NXG=NX;
	NYG=NY;
	NZG=NZ;
	NX = NXG/NPROCX;
	NY = NYG/NPROCY;
	NZ = NZG/NPROCZ;

	nsnap=1+iround((TSNAP2-TSNAP1)/TSNAPINC);

	FP=stdout;

	switch (SNAP) {
		case 1 : /*particle velocity*/

			merge(nsnap,1);
			merge(nsnap,2);
			merge(nsnap,3);
			break;

		case 2 : /*pressure */
			merge(nsnap,6);
			break;

		case 4 : /*particle velocity*/
			merge(nsnap,1);
			merge(nsnap,2);
			merge(nsnap,3);

		case 3 :/*curl and divergence energy*/
			merge(nsnap,4);
			merge(nsnap,5);
			break;

		case 5 :/*Gradient/Model*/
			merge(0,7);
			merge(0,8);
			merge(0,9);
			break;

		default :
			warning(" snapmerge: cannot identify content of snapshot !");
			break;

	}

	return 0;

}
