/// \file isosurfaceD.cxx
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

#include <assert.h>
#include <vector>

#include "isosurfaceD.h"
#include "ijktable.h"

using namespace ISOSURFACE_D;
using namespace IJKTABLE;

/// Generate isosurface in arbitrary dimensions
void ISOSURFACE_D::isosurfaceD
(const int dimension, const ISOSURFACE_TABLE & isotable,
 const SCALAR_TYPE * scalar_grid, const int * grid_length,
 const SCALAR_TYPE isovalue,
 int * & intersected_edge_endpoint, int & num_intersected_edges,
 int * & simplex_vert, int & num_simplices)
  // generate isosurface in arbitrary dimensions
  // dimension = volume dimension
  // isotable = hypercube isosurface table for given dimension
  // scalar_grid[] = array of scalar values
  //   point (x0,x1,x2,...) has scalar value
  //     scalar_grid[x0 + x1*grid_length[0] + 
  //                   x2*grid_length[0]*grid_length[1] + ...]
  // grid_length[i] = grid dimension i
  //   # grid points = 
  //      grid_length[0]*grid_length[1]*grid_length[2]*...
  // isovalue = isosurface scalar value
  // intersected_edge_endpoint = 
  //     array of endpoints of edges intersected by isosurface
  //   allocated and returned by genisosurface
  //   j'th endpoint of edge i = intersected_edge_vert[2*i+j]  (j = 0 or 1)
  //   isosurface vertex i lies on intersected edge i
  // num_intersected_edges = number of intersected edges (returned)
  // simplex_vert = array of isosurface simplex vertices
  //   allocated and returned by genisosurface
  //   # vertices of each simplex = DIMENSION
  //   j'th vertex of simplex i = simplex_vert[DIMENSION*i + j]
  // num_simplices = number of isosurface simplices (returned)
{
  assert(dimension > 0);
  assert(dimension == isotable.Dimension());

  const long num_hypercube_vertices = (1L << dimension);

  assert(isotable.Polyhedron().NumVertices() == num_hypercube_vertices);

  vector<int> elist;
  vector<int> slist;

  // initialize output
  intersected_edge_endpoint = NULL;
  num_intersected_edges = 0;
  simplex_vert = NULL;
  num_simplices = 0;

  long num_grid_vert = get_num_grid_vert(dimension, grid_length);
  if (num_grid_vert == 0)
    return;

  const int maxe = num_grid_vert * dimension;

  // allocate memory for grid edges
  bool * eflag = new bool[maxe];
  int * eindex = new int[maxe];
  for (int ie = 0; ie < maxe; ie++)
    eflag[ie] = false;

  // add intersected edges to elist
  int * coord = new int[dimension];
  for (int iv0 = 0; iv0 < num_grid_vert; iv0++) {
    get_coord(iv0, dimension, grid_length, coord);

    for (int d = 0; d < dimension; d++) {
      if (coord[d]+1 < grid_length[d]) {
	int iv1 = adjacent_grid_vert(iv0, d, grid_length);
	// check if edge(iv0, iv1) is intersected by isosurface
	if ((scalar_grid[iv0] < isovalue && scalar_grid[iv1] >= isovalue) ||
	    (scalar_grid[iv0] >= isovalue && scalar_grid[iv1] < isovalue)) {
	  int ie = get_edge(iv0, iv1, dimension, grid_length);
	  if (eflag[ie] == false) {
	    // add edge (iv0, iv1)
	    int je = elist.size()/2;
	    elist.push_back(iv0);
	    elist.push_back(iv1);
	    eflag[ie] = true;
	    eindex[ie] = je;
	  };
	};
      };
    };
  };

  // compute isosurface
  for (int iv0 = 0; iv0 < num_grid_vert; iv0++) {
    get_coord(iv0, dimension, grid_length, coord);

    // check coord[] + (1,1,...,1) is a grid vertex
    bool flag_vertex = true;
    for (int d = 0; d < dimension; d++) {
      if (coord[d]+1 >= grid_length[d]) {
	flag_vertex = false;
      };
    };
    
    if (flag_vertex) {
      long it = 0;
      for (int j = 0; j < num_hypercube_vertices; j++) {
	int iv1 = get_grid_vert(iv0, j, dimension, grid_length);
	if (scalar_grid[iv1] >= isovalue) {
	  it = (it | (1L << j));
	};
      };

      for (int is = 0; is < isotable.NumSimplices(it); is++) {
	for (int jv = 0; jv < dimension; jv++) {
	  int je = isotable.SimplexVertex(it, is, jv);
	  int jv0 = isotable.Polyhedron().EdgeEndpoint(je, 0);
	  int jv1 = isotable.Polyhedron().EdgeEndpoint(je, 1);
	  int kv0 = get_grid_vert(iv0, jv0, dimension, grid_length);
	  int kv1 = get_grid_vert(iv0, jv1, dimension, grid_length);

	  int ie = get_edge(kv0, kv1, dimension, grid_length);

	  if (eflag[ie] == false) {
	    // programming error
	    // edge ie is not in list of intersected edges
	    delete [] eflag;
	    delete [] eindex;
	    throw ISOSURFACE_TABLE_ERROR("Programming error in isosurface generation. Error computing intersected edges.");
	  };

	  int ke = eindex[ie];
	  slist.push_back(ke);
	};
      };
    };
  };
  delete [] coord;

  delete [] eflag;
  delete [] eindex;

  // allocate and store intersected edges and isosurface simplices
  intersected_edge_endpoint = new int[elist.size()];
  for (int i = 0; i < elist.size(); i++)
    intersected_edge_endpoint[i] = elist[i];
  num_intersected_edges = elist.size()/2;

  simplex_vert = new int[slist.size()];
  for (int i = 0; i < slist.size(); i++)
    simplex_vert[i] = slist[i];
  num_simplices = slist.size()/dimension;
}


/// Return number of grid vertices
long ISOSURFACE_D::get_num_grid_vert
  (const int dimension, const int * grid_length)
  // get number of grid vertices
{
  long num_vert = 1;
  for (int d = 0; d < dimension; d++)
    num_vert = num_vert * grid_length[d];
  return(num_vert);
}

/// Return coordinate of specified vertex
void ISOSURFACE_D::get_coord
  (const int iv, const int dimension, const int * grid_length, int * coord)
  // get coordinates of grid vertex iv
  // Precondition: grid_length[d] > 0 for all d = 0,...,dimension-1
{
  int k = iv;
  for (int d = 0; d < dimension; d++) {
    coord[d] = k % grid_length[d];
    k = k / grid_length[d];
  };
}

/// Return vertex corresponding to specified coordinates
void ISOSURFACE_D::get_vert
  (const int * coord, const int dimension, const int * grid_length, 
   int & iv)
  // get vertex corresponding to coordinates coord[]
{
  iv = 0;
  int k = 1;
  for (int d = 0; d < dimension; d++) {
    iv = iv + coord[d]*k;
    k = k*grid_length[d];
  };
}

int ISOSURFACE_D::get_edge(const int iv0, const int iv1,
			   const int dimension, const int * grid_length)
  // get index of edge (iv0, iv1)
  // Precondition: (iv0, iv1) is a grid edge
{
  int u0, u1;

  if (iv0 < iv1) {
    u0 = iv0;
    u1 = iv1;
  }
  else {
    u0 = iv1;
    u1 = iv0;
  };

  int k0 = u0;
  int k1 = u1;
  for (int d = 0; d < dimension; d++) {
    int c0 = k0 % grid_length[d];
    int c1 = k1 % grid_length[d];
    if (c1 == c0+1) {
      int ie = u0*dimension + d;
      return(ie);
    };
    k0 = k0 / grid_length[d];
    k1 = k1 / grid_length[d];
  };

  // programming error
  throw ISOSURFACE_TABLE_ERROR("Programming error in isosurface generation. Called get_edge for non-grid edge.");
}

int ISOSURFACE_D::adjacent_grid_vert
  (const int iv0, const int d, const int * grid_length)
  // get adjacent vertex in direction d
  // Precondition: iv0 is not the last vertex on the row in direction d
{
  int k = 1;
  for (int i = 0; i < d; i++) {
    k = k * grid_length[i];
  };
  return(iv0+k);
}

int ISOSURFACE_D::get_grid_vert
  (const int iv0, const int j, const int dimension, const int * grid_length)
  // get grid vertex corresponding to j'th hypercube vertex
  // iv0 = grid vertex
  // j = hypercube vertex
{
  int k = 1;
  int j0 = j;
  int iv1 = iv0;
  for (int d = 0; d < dimension; d++) {
    if ((j0 % 2) == 1) {
      iv1 += k;
    };
    j0 = j0/2;
    k = k * grid_length[d];
  };

  return(iv1);
}

