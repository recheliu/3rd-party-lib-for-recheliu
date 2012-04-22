/// \file ijktable.cxx
/// Class containing a table of isosurface patches in a given polyhedron.

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2006, 2003, 2001 Rephael Wenger

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

#include <ctype.h>
#include <limits>
#include <limits.h>
#include <stdlib.h>
// MOD-BY-LEETEN 04/21/2012-FROM:	#include <string.h>
#include <string>
// MOD-BY-LEETEN 04/21/2012-END

#include "ijktable.h"

using namespace IJKTABLE;

#ifndef LONG_BIT

#define LONG_BIT (CHAR_BIT * sizeof(long))

#endif

//**************************************************
// ISOSURFACE_TABLE_POLYHEDRON
//**************************************************

ISOSURFACE_TABLE_POLYHEDRON::ISOSURFACE_TABLE_POLYHEDRON(const int d)
// constructor
{
  Init();
  dimension = d;
}

ISOSURFACE_TABLE_POLYHEDRON::~ISOSURFACE_TABLE_POLYHEDRON()
// destructor
{
  FreeAll();
}

ISOSURFACE_TABLE_POLYHEDRON::ISOSURFACE_TABLE_POLYHEDRON
  (const ISOSURFACE_TABLE_POLYHEDRON & init)
// copy
{
  Init();

  *this = init;
}

const ISOSURFACE_TABLE_POLYHEDRON & ISOSURFACE_TABLE_POLYHEDRON::operator = 
  (const ISOSURFACE_TABLE_POLYHEDRON & right)  
// assign 
{
  if (&right != this) {         // avoid self-assignment
    FreeAll();
    dimension = right.Dimension();
    SetSize(right.NumVertices(), right.NumEdges(), right.NumFacets());

    // copy vertices
    for (int iv = 0; iv < NumVertices(); iv++)
      for (int ic = 0; ic < Dimension(); ic++)
	SetVertexCoord(iv, ic, right.VertexCoord(iv, ic));

    // copy edges
    for (int ie = 0; ie < NumEdges(); ie++)
      SetEdge(ie, right.EdgeEndpoint(ie, 0), right.EdgeEndpoint(ie, 1));

    // copy facets
    for (int jf = 0; jf < NumFacets(); jf++)
      facet[jf] = right.Facet(jf);
  };

  return *this;
}

void ISOSURFACE_TABLE_POLYHEDRON::FreeAll()
// free all memory
{
  num_vertices = 0;
  num_edges = 0;
  num_facets = 0;
  delete [] vertex_coord;
  vertex_coord = NULL;
  delete [] edge_endpoint;
  edge_endpoint = NULL;
  delete [] facet;
  facet = NULL;
}

void ISOSURFACE_TABLE_POLYHEDRON::Init()
// initialize
{
  dimension = 0;
  num_vertices = 0;
  num_edges = 0;
  num_facets = 0;
  vertex_coord = NULL;
  edge_endpoint = NULL;
  facet = NULL;
}

bool ISOSURFACE_TABLE_POLYHEDRON::VertexSign
  (const long code, const int iv) const
// sign of vertex iv in code
{
  long mask = 1L << iv;

  if ((mask & code) == 0)
    return(false);
  else
    return(true);
}

int ISOSURFACE_TABLE_POLYHEDRON::NumFacetVertices(const FACET_INDEX jf) const
// return number of facet vertices
// jf = facet index
{
  int numf = 0;
  for (int iv = 0; iv < NumVertices(); iv++) {
    if (FacetVertexFlag(jf, iv) != 0)
      numf++;
  };

  return(numf);
}

void ISOSURFACE_TABLE_POLYHEDRON::SetDimension(const int d)
// set polyhedron dimension
{
  FreeAll();
  num_vertices = num_edges = 0;

  dimension = d;
}

void ISOSURFACE_TABLE_POLYHEDRON::SetNumVertices(const int numv)
// set number of vertices
// Must be called before setting polyhedron vertices
{
  if (!CheckDimension())
    throw ISOSURFACE_TABLE_ERROR("Illegal polyhedron dimension.");

  FreeAll();
  num_vertices = num_edges = 0;

  if (numv == 0)
    throw ISOSURFACE_TABLE_ERROR("Number of vertices must be non-zero.");

  // Note that even if numv <= LONG_BIT, there may not be enough 
  //   memory to store the isosurface table.
  if (numv > LONG_BIT)
    throw ISOSURFACE_TABLE_ERROR
      ("Number of polyhedron vertices is too large.");

  num_vertices = numv;
  vertex_coord = new int[num_vertices*Dimension()];
}

void ISOSURFACE_TABLE_POLYHEDRON::SetNumEdges(const int nume)
// set number of edges
// Must be called before setting polyhedron edges
{
  delete [] edge_endpoint;
  edge_endpoint = NULL;
  num_edges = 0;

  if (!CheckDimension())
    throw ISOSURFACE_TABLE_ERROR("Illegal dimension.");

  if (NumVertices() == 0)
    throw ISOSURFACE_TABLE_ERROR
      ("Number of vertices must be set before number of edges.");

  if (nume < 1)
    throw ISOSURFACE_TABLE_ERROR("Number of edges must be non-zero.");

  if (nume > std::numeric_limits<EDGE_INDEX>::max())
    throw ISOSURFACE_TABLE_ERROR("Number of polyhedron edges is too large.");

  num_edges = nume;
  edge_endpoint = new int[num_edges*2];
}


void ISOSURFACE_TABLE_POLYHEDRON::SetNumFacets(const int numf)
// set number of facets
// Must be called before setting polyhedron facets
{
  delete [] facet;
  facet = NULL;
  num_facets = 0;

  if (!CheckDimension())
    throw ISOSURFACE_TABLE_ERROR("Illegal dimension.");

  if (NumVertices() == 0)
    throw ISOSURFACE_TABLE_ERROR
      ("Number of vertices must be set before number of facets.");

  if (numf < 1)
    throw ISOSURFACE_TABLE_ERROR("Number of facets must be non-zero.");

  if (numf > std::numeric_limits<FACET_INDEX>::max())
    throw ISOSURFACE_TABLE_ERROR("Number of polyhedron facets is too large.");

  num_facets = numf;
  facet = new FACET[numf];

  if (facet == NULL)
    throw ISOSURFACE_TABLE_ERROR
      ("Unable to allocate memory for list of facets.");

  // initialize each facet to 0
  for (int jf = 0; jf < numf; jf++) {
    facet[jf] = 0;
  };

}

void ISOSURFACE_TABLE_POLYHEDRON::SetVertexCoord
  (const int iv, const int ic, const int coord)
// set polyhedron vertex coordinate
// iv = vertex index.  In range [0..NumVertices()-1].
// ic = coordinate index. In range [0..Dimension()-1].
// coord = coordinate.  Must be even.
{
  if (iv < 0 || iv >= NumVertices())
    throw ISOSURFACE_TABLE_ERROR("Illegal vertex index.");

  if (ic < 0 || ic >= Dimension())
    throw ISOSURFACE_TABLE_ERROR("Illegal vertex coordinate index.");

  if (coord%2 != 0)
    throw ISOSURFACE_TABLE_ERROR
      ("Illegal vertex coordinate.  Vertex coordinate must be even.");

  if (vertex_coord == NULL)
    throw ISOSURFACE_TABLE_ERROR
      ("Vertex coordinate memory not allocated.");

  vertex_coord[iv*Dimension() + ic] = coord;
}

void ISOSURFACE_TABLE_POLYHEDRON::SetEdge
  (const EDGE_INDEX ie, const int iv0, const int iv1)
// set polyhedron edge coordinate
// ie = edge index.  In range [0..NumEdges()-1].
// iv0 = vertex 0 index.  In range [0..NumVertices()-1].
// iv1 = vertex 1 index.  In range [0..NumVertices()-1].
{
  int ie_int = int(ie);
  if (ie_int < 0 || ie_int >= NumEdges())
    throw ISOSURFACE_TABLE_ERROR("Illegal edge index.");

  if (iv0 < 0 || iv0 > NumVertices() ||
      iv1 < 0 || iv1 > NumVertices())
    throw ISOSURFACE_TABLE_ERROR("Illegal vertex index.");

  if (edge_endpoint == NULL)
    throw ISOSURFACE_TABLE_ERROR("Edge endpoint memory not allocated.");

  edge_endpoint[2*int(ie)] = iv0;
  edge_endpoint[2*int(ie)+1] = iv1;
}

void ISOSURFACE_TABLE_POLYHEDRON::SetFacetVertex
  (const FACET_INDEX jf, const int iv, const bool in_facet)
// set facet vertex to in or not in facet
// jf = facet index
// iv = vertex index
// in_facet = true if vertex iv is in facet
// Note: SetNumFacets must be called before SetFacetVertex
{
  int jf_int = int(jf);
  if (jf_int < 0 || jf_int >= NumFacets())
    throw ISOSURFACE_TABLE_ERROR("Illegal facet index.");

  if (iv < 0 || iv >= NumVertices()) {
    throw ISOSURFACE_TABLE_ERROR("Illegal vertex index.");
  };

  if (facet == NULL)
    throw ISOSURFACE_TABLE_ERROR
      ("Number of facets must be set before facet vertices.");

  long mask = 1L << iv;
  if (in_facet) {
    facet[jf_int] = facet[jf_int] | mask;
  }
  else {
    mask = ~mask;
    facet[jf_int] = facet[jf_int] & mask;
  };
}

int ISOSURFACE_TABLE_POLYHEDRON::MidpointCoord
  (const EDGE_INDEX ie, const int ic) const
// return ic'th coordinate of midpoint of edge ie
// Note: vertice coordinates must all be even so midpoint coordinate
//   is an integer
{
  int iv0 = EdgeEndpoint(ie, 0);
  int iv1 = EdgeEndpoint(ie, 1);
  int coord0 = VertexCoord(iv0, ic);
  int coord1 = VertexCoord(iv1, ic);
  
  return((coord0+coord1)/2);
}

void ISOSURFACE_TABLE_POLYHEDRON::GenCube(const int cube_dimension)
// generate a polyhedron
{
  dimension = cube_dimension;

  if (!CheckDimension())
    throw ISOSURFACE_TABLE_ERROR();

  int numv = 1L << Dimension();
  int nume = (numv * Dimension()) / 2;
  int numf = 2*Dimension();

  SetSize(numv, nume, numf);

  // set vertex coordinates
  for (int iv = 0; iv < NumVertices(); iv++) {
    long mask = 1L;
    for (int ic = 0; ic < Dimension(); ic++) {
      int bit = mask & iv;
      int coord = 0;
      if (bit != 0)
	coord = 2;
      SetVertexCoord(iv, ic, coord);
      mask = mask << 1;
    };
  };

  // generate edges in lexicographic order
  int ie = 0;
  long control = 0;
  while (ie < NumEdges()) {

    // find first 0 bit in control
    int ic = 0;
    long mask = 1L;
    while ((mask & control) != 0) {
      mask = mask << 1;
      ic++;
    };

    // find start vertex by stripping ic bits from control
    int start = control;
    start = start >> ic;
    start = start << ic;

    for (int i = 0; i < (1L << ic); i++) {
      int iv0 = start + i;
      int iv1 = iv0 + (1L << ic);
      SetEdge(ie, iv0, iv1);
      ie++;
    };

    control++;
  };

  if (control+1 != (1L << Dimension()))
    throw ISOSURFACE_TABLE_ERROR("Programming error detected in GenCube.");

  // generate facets
  int num_vertices_per_facet = NumVertices()/2;
  long mask = 1L;
  for (int ic = 0; ic < Dimension(); ic++) {
    int jf0 = 2*ic;
    int jf1 = 2*ic+1;

    for (int iv = 0; iv < NumVertices(); iv++) {
      long bit = mask & iv;
      if (bit == 0) {
	SetFacetVertex(jf0, iv, true);
      }
      else {
	SetFacetVertex(jf1, iv, true);
      };
    };

    if (NumFacetVertices(jf0) != num_vertices_per_facet ||
	NumFacetVertices(jf1) != num_vertices_per_facet) {
      throw ISOSURFACE_TABLE_ERROR("Programming error detected in GenCube.");
    };

    mask = mask << 1;
  };

}


void ISOSURFACE_TABLE_POLYHEDRON::GenSimplex(const int simplex_dimension)
// generate a simplex
{
  dimension = simplex_dimension;

  if (!CheckDimension())
    throw ISOSURFACE_TABLE_ERROR();

  int numv = Dimension() + 1;
  int nume = (numv * Dimension()) / 2;
  int numf = Dimension() + 1;

  SetSize(numv, nume, numf);

  // initialize all vertex coordinates to 0
  for (int iv = 0; iv < NumVertices(); iv++)
    for (int ic = 0; ic < Dimension(); ic++)
      SetVertexCoord(iv, ic, 0);

  // set vertex coordinates
  int ic = 0;
  for (int jv = 1; jv < NumVertices(); jv++) {
    SetVertexCoord(jv, ic, 2);
    ic++;
  };

  // generate edges in lexicographic order
  int ie = 0;
  for (int iv0 = 0; iv0 < NumVertices()-1; iv0++)
    for (int iv1 = iv0+1; iv1 < NumVertices(); iv1++) {
      SetEdge(ie, iv0, iv1);
      ie++;
    };

  if (ie != nume)
    throw ISOSURFACE_TABLE_ERROR("Programming error detected in GenSimplex.");

  // generate facets
  int num_vertices_per_facet = Dimension();
  for (int jf = 0; jf < numf; jf++) {

    for (int jv = 0; jv < num_vertices_per_facet; jv++) {
      int iv = (jf + jv) % numv;
      SetFacetVertex(jf, iv, true);
    };

    if (NumFacetVertices(jf) != num_vertices_per_facet) 
      throw ISOSURFACE_TABLE_ERROR
	("Programming error detected in GenSimplex.");
  };
}


void ISOSURFACE_TABLE_POLYHEDRON::WriteVertexCoordText(ostream & out) const
// write vertex coordinates in text format
{
  for (int iv = 0; iv < NumVertices(); iv++) {
    for (int ic = 0; ic < Dimension(); ic++) {
      out << "  " << VertexCoord(iv, ic);
    };
    out << endl;
  };
}

void ISOSURFACE_TABLE_POLYHEDRON::WriteEdgeListText(ostream & out) const
// write edge list in text format
{
  for (int ie = 0; ie < NumEdges(); ie++) {
    int iv0 = EdgeEndpoint(ie, 0);
    int iv1 = EdgeEndpoint(ie, 1);
    out << "  " << iv0 << "  " << iv1 << endl;
  };
}

void ISOSURFACE_TABLE_POLYHEDRON::WriteFacetVerticesText(ostream & out) const
// write facet vertices lists in text format
{
  for (int jf = 0; jf < NumFacets(); jf++) {
    int numv = NumFacetVertices(jf);
    out << "  " << numv << endl;
    int numv2 = 0;
    for (int iv = 0; iv < NumVertices(); iv++) {
      if (FacetVertexFlag(jf, iv)) {
	out << "  " << iv;
	numv2++;
      };
    };
    if (numv != numv2) {
      throw ISOSURFACE_TABLE_ERROR
	("Error in number of facet vertices.");
    };
    out << endl;
  };
}

void ISOSURFACE_TABLE_POLYHEDRON::ReadVertexCoordText(istream & in)
// read vertex coordinates in text format
{
  for (int iv = 0; iv < NumVertices(); iv++) {
    for (int ic = 0; ic < Dimension(); ic++) {
      if (in.eof())
	throw ISOSURFACE_TABLE_ERROR
	  ("Premature end of file while reading polyhedron vertices.");
      if (!in)
	throw ISOSURFACE_TABLE_ERROR("Error reading polyhedron vertices.");

      int coord;
      in >> coord;

      if (coord%2 != 0)
	throw ISOSURFACE_TABLE_ERROR
	  ("Polyhedron coordinates must be even.");

      SetVertexCoord(iv, ic, coord);
    };
  };
}

void ISOSURFACE_TABLE_POLYHEDRON::ReadEdgeListText(istream & in)
// read edge endpoints in text format
{
  for (int ie = 0; ie < NumEdges(); ie++) {
    if (in.eof())
      throw ISOSURFACE_TABLE_ERROR
	("Premature end of file while reading polyhedron edges.");
    if (!in)
      throw ISOSURFACE_TABLE_ERROR("Error reading polyhedron edges.");

    int iv0, iv1;
    in >> iv0 >> iv1;

    if (iv0 < 0 || iv1 < 0 || iv0 >= NumVertices() || iv1 >= NumVertices())
      throw ISOSURFACE_TABLE_ERROR
	("Illegal vertex index in polyhedron edge list.");

    SetEdge(ie, iv0, iv1);
  };
}

void ISOSURFACE_TABLE_POLYHEDRON::ReadFacetsText(istream & in)
// read lists of facet vertices text format
{
  for (int jf = 0; jf < NumFacets(); jf++) {
    if (in.eof())
      throw ISOSURFACE_TABLE_ERROR
	("Premature end of file while reading polyhedron facets.");
    if (!in)
      throw ISOSURFACE_TABLE_ERROR("Error reading polyhedron facets.");

    int num_facet_vertices;
    in >> num_facet_vertices;
    if (num_facet_vertices < 0 || num_facet_vertices > NumVertices()) {
      throw ISOSURFACE_TABLE_ERROR
	("Illegal number of polyhedron facet vertices.");
    };

    for (int jv = 0; jv < num_facet_vertices; jv++) {
      int iv;
      in >> iv;
      if (iv < 0 || iv >= NumVertices()) {
	throw ISOSURFACE_TABLE_ERROR
	  ("Illegal vertex index in polyhedron facet.");
      };

      SetFacetVertex(jf, iv, true);
    };
  };
}

bool ISOSURFACE_TABLE_POLYHEDRON::CheckDimension() const
// check dimension
{
  if (dimension < 1)
    return(false);
  else
    return(true);
}

bool ISOSURFACE_TABLE_POLYHEDRON::Check() const
// check polyhedron
{
  if (!CheckDimension())
    return(false);

  if (NumVertices() < 1 || NumEdges() < 1 ||
      vertex_coord == NULL ||
      edge_endpoint == NULL)
    return(false);

  for (int iv = 0; iv < NumVertices(); iv++) {
    for (int ic = 0; ic < Dimension(); ic++) {
      if ((VertexCoord(iv, ic) % 2) != 0)
	return(false);
    };
  };

  for (int ie = 0; ie < NumEdges(); ie++) {
    for (int ip = 0; ip < 2; ip++) {
      int iv = EdgeEndpoint(ie, ip);
      if (iv < 0 || iv >= NumVertices())
	return(false);
    };
  };

  if (NumFacets() > 0) {
    if (facet == NULL)
      return(false);
  };

  return(true);
}

//**************************************************
// ISOSURFACE_TABLE
//**************************************************

ISOSURFACE_TABLE::ISOSURFACE_TABLE_ENTRY::ISOSURFACE_TABLE_ENTRY()
// constructor
{
  num_simplices = 0;
  simplex_vertex_list = NULL;
}

ISOSURFACE_TABLE::ISOSURFACE_TABLE_ENTRY::~ISOSURFACE_TABLE_ENTRY()
// destructor
{
  delete [] simplex_vertex_list;
  simplex_vertex_list = NULL;
  num_simplices = 0;
}

bool ISOSURFACE_TABLE::ISOSURFACE_TABLE_ENTRY::Check() const
{
  if (num_simplices < 0)
    return(false);

  if (num_simplices > 0 && simplex_vertex_list == NULL)
    return(false);

  return(true);
}

// initializer for isosurface table file header
const char * const ISOSURFACE_TABLE::FILE_HEADER = "<ISOTABLE2>";
const char * const ISOSURFACE_TABLE::VERTEX_HEADER = "<VERTICES>";
const char * const ISOSURFACE_TABLE::EDGE_HEADER = "<EDGES>";
const char * const ISOSURFACE_TABLE::FACET_HEADER = "<FACETS>";
const char * const ISOSURFACE_TABLE::TABLE_HEADER = "<TABLE>";
const char * const ISOSURFACE_TABLE::FILE_HEADER_V1 = "<ISOTABLE>";

ISOSURFACE_TABLE::ISOSURFACE_TABLE(const int d) :
  polyhedron(d)
// constructor
// d = dimension of space containing isosurface.  Should be 2, 3 or 4.
{
  max_num_vertices = 20;
  // Note: Even tables for polyhedron of this size are probably impossible 
  //   to compute/store

  num_table_entries = 0;
  entry = NULL;
  is_table_allocated = false;
  if (!CheckDimension())
    throw ISOSURFACE_TABLE_ERROR();
}

ISOSURFACE_TABLE::~ISOSURFACE_TABLE()
// destructor
{
  num_table_entries = 0;
  delete [] entry;
  entry = NULL;
}

void ISOSURFACE_TABLE::AllocTable()
// allocate memory for isosurface table
{
  if (Polyhedron().Dimension() == 0 ||
      Polyhedron().NumVertices() == 0 ||
      Polyhedron(). NumEdges() == 0) {
    throw ISOSURFACE_TABLE_ERROR
      ("Create polyhedron before allocating isosurface table.");
  };
      
  if (IsTableAllocated() || entry != NULL) {
    throw ISOSURFACE_TABLE_ERROR
      ("Isosurface table cannot be allocated more than once.");
  };

  delete [] entry;
  entry = NULL;
  num_table_entries = 0;

  if (Polyhedron().NumVertices() > MaxNumVertices())
    throw ISOSURFACE_TABLE_ERROR
      ("Table polyhedron contains too many vertices.");

  // use shift operator to calculate 2^(Polyhedron().NumVertices())
  num_table_entries = 1L << Polyhedron().NumVertices();
  entry = new ISOSURFACE_TABLE_ENTRY[num_table_entries];

  is_table_allocated = true;
}


void ISOSURFACE_TABLE::CheckTableAlloc() const
// throw error if table is allocated
{
  if (IsTableAllocated()) {
    throw ISOSURFACE_TABLE_ERROR
      ("Illegal attempt to define polyhedron after allocating isosurface table.");
  };
}

void ISOSURFACE_TABLE::SetNumSimplices(const TABLE_INDEX it, const int nums)
// set number of simplices in table entry it
// it = table entry
// nums = number of simplices
{
  if (!IsTableAllocated() || entry == NULL) {
    throw ISOSURFACE_TABLE_ERROR
      ("Table must be allocated before entering table entries.");
  };

  if (it < 0 || it >= NumTableEntries() ||
      nums < 0)
    throw ISOSURFACE_TABLE_ERROR();

  entry[it].num_simplices = 0;
  delete entry[it].simplex_vertex_list;
  entry[it].simplex_vertex_list = NULL;

  if (nums > 0)
    entry[it].simplex_vertex_list = 
      new EDGE_INDEX[nums*NumVerticesPerSimplex()];

  entry[it].num_simplices = nums;
}


void ISOSURFACE_TABLE::SetSimplexVertex
  (const TABLE_INDEX it, const int is, const int iv, const EDGE_INDEX ie)
// set simplex vertex
// it = index table entry.  In range [0..NumTableEntries()-1].
// is = index simplex.  
// iv = index simplex vertex.  In range [0..NumVerticesPerSimplex()-1].
// ie = index of edge containing simplex vertex
{
  entry[it].simplex_vertex_list[is*NumVerticesPerSimplex()+iv] = ie;
}

void ISOSURFACE_TABLE::WriteText(ostream & out)
// write polyhedron & table, text format
{
  WriteFileHeaderText(out);
  out << Dimension() << endl;
  out << VertexHeader() << endl;
  out << polyhedron.NumVertices() << endl;
  polyhedron.WriteVertexCoordText(out);
  out << EdgeHeader() << endl;
  out << polyhedron.NumEdges() << endl;
  polyhedron.WriteEdgeListText(out);
  if (polyhedron.NumFacets() > 0) {
    out << FacetHeader() << endl;
    out << polyhedron.NumFacets() << endl;
    polyhedron.WriteFacetVerticesText(out);
  };
  out << TableHeader() << endl;
  WriteTableText(out);
}

void ISOSURFACE_TABLE::WriteTextV1(ostream & out)
// write polyhedron & table, text format
// old version
{
  out << FileHeaderV1() << endl;
  out << Dimension() << endl;
  out << polyhedron.NumVertices() << endl;
  polyhedron.WriteVertexCoordText(out);
  out << polyhedron.NumEdges() << endl;
  polyhedron.WriteEdgeListText(out);
  WriteTableText(out);
}

void ISOSURFACE_TABLE::WriteFileHeaderText(ostream & out)
// write isosurface table file header
{
  out << FileHeader() << endl;
}

void ISOSURFACE_TABLE::WriteTableText(ostream & out)
// write isosurface table in text format
{
  for (int it = 0; it < NumTableEntries(); it++) {
    out << NumSimplices(it) << endl;
    for (int is = 0; is < NumSimplices(it); is++) {
      for (int iv = 0; iv < NumVerticesPerSimplex(); iv++) {
	out << "  " << int(SimplexVertex(it, is, iv));
      };
      out << endl;
    };
  };
}

void ISOSURFACE_TABLE::ReadText(istream & in)
// read polyhedron & table, text format
{
  int d;
  int num_poly_vert;
  int num_poly_edges;
  int num_poly_facets;
  int num_facet_vertices;
  std::string s;

  if (!in)
    throw ISOSURFACE_TABLE_ERROR("Input stream not open.");

  in >> s;
  if (s.compare(FILE_HEADER_V1) == 0) {
    ReadNoHeaderTextV1(in);
    return;
  }
  else if (s.compare(FILE_HEADER) != 0) {
    throw ISOSURFACE_TABLE_ERROR("Illegal file header format.");
  };

  in >> d;
  if (d < 2)
    throw ISOSURFACE_TABLE_ERROR("Dimension must be at least 2.");

  polyhedron.SetDimension(d);

  in >> s;
  if (s.compare(VERTEX_HEADER) != 0) {
    throw ISOSURFACE_TABLE_ERROR("Illegal file vertex format.");
  };

  in >> num_poly_vert;
  if (num_poly_vert < 1)
    throw ISOSURFACE_TABLE_ERROR
      ("Number of polyhedron vertices must be at least 1.");
    
  polyhedron.SetNumVertices(num_poly_vert);
  polyhedron.ReadVertexCoordText(in);

  in >> s;
  if (s.compare(EDGE_HEADER) != 0) {
    throw ISOSURFACE_TABLE_ERROR("Illegal file edge format.");
  };

  in >> num_poly_edges;
  if (num_poly_edges < 1)
    throw ISOSURFACE_TABLE_ERROR
      ("Number of polyhedron edges must be at least 1.");

  polyhedron.SetNumEdges(num_poly_edges);
  polyhedron.ReadEdgeListText(in);

  in >> s;
  if (s.compare(FACET_HEADER) == 0) {

    in >> num_poly_facets;
    if (num_poly_facets < 1)
      throw ISOSURFACE_TABLE_ERROR
	("Number of polyhedron facets must be at least 1.");

    polyhedron.SetNumFacets(num_poly_facets);
    polyhedron.ReadFacetsText(in);

    // read next header
    in >> s;
  };

  AllocTable();

  if (s.compare(TABLE_HEADER) != 0) {
    throw ISOSURFACE_TABLE_ERROR("Illegal file table format.");
  };

  ReadTableText(in);
}

void ISOSURFACE_TABLE::ReadTextV1(istream & in)
// read polyhedron & table, old version text format
{
  std::string s;

  if (!in)
    throw ISOSURFACE_TABLE_ERROR("Input stream not open.");

  in >> s;
  if (s.compare(FILE_HEADER_V1) != 0) {
    throw ISOSURFACE_TABLE_ERROR("Illegal file header format.");
  };

  ReadNoHeaderTextV1(in);
}


void ISOSURFACE_TABLE::ReadNoHeaderTextV1(istream & in)
// read polyhedron & table, old version text format
// do not read file header
{
  int d;
  int num_poly_vert;
  int num_poly_edges;
  std::string s;

  if (!in)
    throw ISOSURFACE_TABLE_ERROR("Input stream not open.");

  in >> d;
  if (d < 2)
    throw ISOSURFACE_TABLE_ERROR("Dimension must be at least 2.");

  polyhedron.SetDimension(d);

  in >> num_poly_vert;
  if (num_poly_vert < 1)
    throw ISOSURFACE_TABLE_ERROR
      ("Number of polyhedron vertices must be at least 1.");
    
  polyhedron.SetNumVertices(num_poly_vert);
  polyhedron.ReadVertexCoordText(in);

  in >> num_poly_edges;
  if (num_poly_edges < 1)
    throw ISOSURFACE_TABLE_ERROR
      ("Number of polyhedron edges must be at least 1.");

  polyhedron.SetNumEdges(num_poly_edges);
  polyhedron.ReadEdgeListText(in);

  AllocTable();

  ReadTableText(in);
}


void ISOSURFACE_TABLE::ReadFileHeaderText(istream & in)
// read isosurface table file header in text format
{
  char c;

  if (!in)
    throw ISOSURFACE_TABLE_ERROR
      ("Error reading isosurface table.");

  while (in.get(c))
    {
      if (!isspace(c)) {
	break;
      };
    };

  if (in.eof())
    throw ISOSURFACE_TABLE_ERROR("Empty file.");
  if (!in)
    throw ISOSURFACE_TABLE_ERROR
      ("Error reading isosurface table.");

  int buffer_length = strlen(FileHeader())+1;
  char * buffer = new char[buffer_length];

  int i = 0;
  while ((i < buffer_length-1) && !isspace(c) && in) {
    buffer[i] = c;
    i++;
    if (!in.get(c))
      break;
  };
  buffer[i] = '\0';

  if (strcmp(FileHeader(), buffer) != 0) {
    delete [] buffer;
    throw ISOSURFACE_TABLE_ERROR("Illegal file header.");
  };

  delete [] buffer;
}

void ISOSURFACE_TABLE::ReadTableText(istream & in)
// read isosurface table in text format
{
  for (int it = 0; it < NumTableEntries(); it++) {

    if (in.eof())
      throw ISOSURFACE_TABLE_ERROR("Premature end of file.");
    if (!in)
      throw ISOSURFACE_TABLE_ERROR
	("Error reading isosurface table.");

    int nums;
    in >> nums;
    if (nums < 0)
      throw ISOSURFACE_TABLE_ERROR
	("Number of simplices cannot be negative.");

    SetNumSimplices(it, nums);
    for (int is = 0; is < NumSimplices(it); is++) {
      for (int iv = 0; iv < NumVerticesPerSimplex(); iv++) {
	int ie;
	in >> ie;
	if (ie < 0 || ie >= polyhedron.NumEdges())
	  throw ISOSURFACE_TABLE_ERROR
	    ("Illegal edge index in isosurface table.");

	SetSimplexVertex(it, is, iv, EDGE_INDEX(ie));
      };
    };
  };
}

bool ISOSURFACE_TABLE::CheckDimension(const int d) const
// check dimension
{
  if (d < 1)
    return(false);
  else
    return(true);
}

bool ISOSURFACE_TABLE::CheckTable() const
// check table
{
  if (polyhedron.NumVertices() >= LONG_BIT || polyhedron.NumVertices() < 1)
    return(false);

  if (polyhedron.NumVertices() > MaxNumVertices())
    return(false);

  if (NumTableEntries() != (1L << polyhedron.NumVertices()))
    return(false);

  if (entry == NULL)
    return(false);

  for (int it = 0; it < NumTableEntries(); it++)
    if (!entry[it].Check())
      return(false);

  for (int jt = 0; jt < NumTableEntries(); jt++)
    for (int is = 0; is < NumSimplices(jt); is++)
      for (int iv = 0; iv < NumVerticesPerSimplex(); iv++) {
	int ie = int(SimplexVertex(jt, is, iv));
	if (ie < 0 || ie >= polyhedron.NumEdges())
	  return(false);
      };

  return(true);
}

//**************************************************
// Routines for generating polyhedra
//**************************************************

void IJKTABLE::generate_prism
  (const ISOSURFACE_TABLE_POLYHEDRON & base_polyhedron,
   ISOSURFACE_TABLE_POLYHEDRON & prism)
// generate a prism with base base_polyhedron
// first numv vertices have last coordinate = 0
// last numv vertices have last coordinate = 2
// first nume edges connect first numv vertices
// second nume edges connect second numv vertices
// third numv edges connect first to second set of vertices
// (numv = # vertices in base_polyhedron; nume = # edges in base_polyhedron)
{
  int dim = base_polyhedron.Dimension();
  int numc = dim;
  int numv = base_polyhedron.NumVertices();
  int nume = base_polyhedron.NumEdges();
  int numf = base_polyhedron.NumFacets();
  int prism_dim = dim + 1;
  int prism_numc = prism_dim;
  int prism_lastc = prism_numc - 1;
  int prism_numv = numv * 2;
  int prism_nume = nume * 2 + numv;
  int prism_numf = 2 + numf;
  prism.SetDimension(prism_dim);
  prism.SetSize(prism_numv, prism_nume, prism_numf);

  // set prism vertex coord
  for (int iv = 0; iv < numv; iv++) {
    for (int ic = 0; ic < prism_lastc; ic++) {
      int coord = base_polyhedron.VertexCoord(iv, ic);
      prism.SetVertexCoord(iv, ic, coord);
      prism.SetVertexCoord(iv+numv, ic, coord);
    };
    prism.SetVertexCoord(iv, prism_lastc, 0);
    prism.SetVertexCoord(iv+numv, prism_lastc, 2);
  };

  // set edges
  for (int ie = 0; ie < base_polyhedron.NumEdges(); ie++) {
    int iv0 = base_polyhedron.EdgeEndpoint(ie, 0);
    int iv1 = base_polyhedron.EdgeEndpoint(ie, 1);
    prism.SetEdge(ie, iv0, iv1);
    prism.SetEdge(ie+nume, iv0+numv, iv1+numv);
  };

  for (int iv = 0; iv < base_polyhedron.NumVertices(); iv++) {
    int ie = 2*nume + iv;
    prism.SetEdge(ie, iv, iv+numv);
  };

  // set facets
  for (int iv = 0; iv < numv; iv++) {
    // set two base facets
    prism.SetFacetVertex(0, iv, true);
    prism.SetFacetVertex(1, iv+numv, true);
  };

  for (int jf = 0; jf < numf; jf++) {
    // set prism facet corresponding to original facet jf
    int prism_jf = 2 + jf;
    for (int iv = 0; iv < numv; iv++) {
      if (base_polyhedron.FacetVertexFlag(jf, iv)) {
	prism.SetFacetVertex(prism_jf, iv, true);
	prism.SetFacetVertex(prism_jf, iv+numv, true);
      };
    };
  };

}


