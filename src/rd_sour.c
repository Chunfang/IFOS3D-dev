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
 *   Read external source wavelet 
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float *rd_sour(int *nts,FILE* fp_source){

	/* local variables */
	float *psource;
	int i, c;
	extern int NT;

	if (fp_source==NULL) err(" Source file could no be opened !");
	/* fscanf(fp_source,"%i", nts); */
        *nts=0;
        while ((c=fgetc(fp_source)) != EOF)
         if (c=='\n') ++(*nts);
        rewind(fp_source);
	psource=vector(1,NT);
	for (i=1;i<=*nts;i++) fscanf(fp_source,"%e",&psource[i]);
	fclose(fp_source);
	return psource;
}
