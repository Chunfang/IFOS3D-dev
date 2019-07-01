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

/* -----------------------------------------------------------------------
 * Initialise gradient with zero
 -------------------------------------------------------------------------*/

#include "fd.h"

void zero_grad(int NX, int NY, int NZ, float *** grad1, float *** grad2, float *** grad3){

   int i, j, k;

    for (j=1;j<=NY;j++){
      for (i=1;i<=NX;i++){
	for (k=1;k<=NZ;k++){
	  grad1[j][i][k]=0.0;
	  grad2[j][i][k]=0.0;
	  grad3[j][i][k]=0.0; 
	}
      }
    }

}