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
 *   Write program name etc to stdout                          
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/* printing all important parameters to FILE *fp */
void info(FILE *fp){

	fprintf(fp," \n\n\n***********************************************************\n");
	fprintf(fp," This is program IFOS3D.          \n");
	fprintf(fp," 3D-Adjoint Full WAveform Inversion Tool      \n");
	fprintf(fp," for 3D elastic FWI                                                         \n");
	fprintf(fp," FWI code written by  S. Butzer        \n");
	fprintf(fp," forward FD code  written by  T. Bohlen    \n");
	fprintf(fp," Geophysical Institute, Department of Physics,         \n");
	fprintf(fp," Institute of Technology, Karlsruhe, Germany         \n");
	fprintf(fp," http://www.gpi.kit.edu \n");
	fprintf(fp," ***********************************************************\n");
	fprintf(fp,"\n");
}
