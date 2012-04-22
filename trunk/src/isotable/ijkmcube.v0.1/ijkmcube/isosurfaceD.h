/// \file isosurfaceD.h
/// Generate isosurface in arbitrary dimensions

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2006 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _ISOSURFACE_D_
#define _ISOSURFACE_D_

#include "ijktable.h"

namespace ISOSURFACE_D {

typedef float SCALAR_TYPE;

long get_num_grid_vert(const int dimension, const int * grid_length);
void get_coord(const int iv, const int dimension, const int * grid_length, 
	       int * coord);
void get_vert(const int * coord, const int dimension, const int * grid_length, 
	      int & iv);
int get_edge(const int iv0, const int iv1, 
	     const int dimension, const int * grid_length);
int adjacent_grid_vert(const int iv0, const int d, const int * grid_length);
int get_grid_vert(const int iv0, const int j, 
		  const int dimension, const int * grid_length);

void isosurfaceD
(const int dimension, const IJKTABLE::ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE * scalar_grid, const int * grid_length,
 const SCALAR_TYPE isovalue,
 int * & intersected_edge_endpoint, int & num_intersected_edges,
 int * & simplex_vert, int & num_simplices);

};

#endif
