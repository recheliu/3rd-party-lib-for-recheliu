/// \file ijktable.h
/// Class containing a table of isosurface patches in a given polyhedron.
/// All 2^numv +/- patterns are stored in the table 
///   where numv = # polyhedron vertices.

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

#ifndef _IJKTABLE_
#define _IJKTABLE_

#include <iostream>

using namespace std;

namespace IJKTABLE {

typedef unsigned char u_char;	// ADD-BY-LEETEN 04/21/2012

typedef u_char EDGE_INDEX;
typedef u_char FACET_INDEX;
typedef int TABLE_INDEX;
typedef int FACET;

const int NO_VERTEX = -1;

class ISOSURFACE_TABLE_POLYHEDRON {

 protected:
  int dimension;
  int num_vertices;
  int num_edges;
  int num_facets;
  int * vertex_coord;
  int * edge_endpoint;
  FACET * facet;
  void FreeAll();        // free all memory
  void Init();           // initialize

 public:
  ISOSURFACE_TABLE_POLYHEDRON(const int d);  // constructor
  ~ISOSURFACE_TABLE_POLYHEDRON();            // destructor
  ISOSURFACE_TABLE_POLYHEDRON(const ISOSURFACE_TABLE_POLYHEDRON & init);
  // copy
  const ISOSURFACE_TABLE_POLYHEDRON & operator = 
    (const ISOSURFACE_TABLE_POLYHEDRON &);  // assign 


  // get functions
  int Dimension() const { return(dimension); };
  int NumVertices() const { return(num_vertices); };
  int NumEdges() const { return(num_edges); };
  int NumFacets() const { return(num_facets); };
  int NumFacetVertices(const FACET_INDEX jf) const;
    // jf = facet index
  int VertexCoord(const int iv, const int ic) const
    // iv = vertex index. ic = coordinate index.
    { return(vertex_coord[iv*dimension + ic]); };
  int EdgeEndpoint(const EDGE_INDEX ie, const int j) const
    // ie = edge index. j = 0 or 1.
    { return(edge_endpoint[int(ie)*2 + j]); };
  FACET Facet(const FACET_INDEX jf) const
    // jf = facet index
    { return(facet[jf]); };
  int MidpointCoord(const EDGE_INDEX ie, const int ic) const;
  bool VertexSign(const long code, const int iv) const;
  bool FacetVertexFlag(const FACET_INDEX jf, const int iv) const
    { return(facet[jf] & ((1L) << iv)); };

  // set functions
  void SetDimension(const int d);
  void SetNumVertices(const int numv);
  void SetNumEdges(const int nume);
  void SetNumFacets(const int numf);
  void SetSize(const int numv, const int nume)
    { SetNumVertices(numv); SetNumEdges(nume); };
  void SetSize(const int numv, const int nume, const int numf)
    { SetNumVertices(numv); SetNumEdges(nume); SetNumFacets(numf); };
  void SetVertexCoord(const int iv, const int ic, const int coord);
  // Note: SetNumVertices or SetSize must be called before SetVertexCoord
  void SetEdge(const EDGE_INDEX ie, const int iv0, const int iv1);
  // Note: SetNumEdges or SetSize must be called before SetEdge
  void SetFacetVertex(const FACET_INDEX jf, const int iv, const bool in_facet);
  // Note: SetNumFacets must be called before SetFacetVertex

  // generate polyhedron
  void GenCube(const int cube_dimension);
  void GenSimplex(const int simplex_dimension);

  // write polyhedron information
  void WriteVertexCoordText(ostream & out) const;
  void WriteEdgeListText(ostream & out) const;
  void WriteFacetVerticesText(ostream & out) const;

  // read polyhedron information
  void ReadVertexCoordText(istream & in);
  void ReadEdgeListText(istream & in);
  void ReadFacetsText(istream & in);

  // check functions
  bool CheckDimension() const;
  bool Check() const;
};

typedef ISOSURFACE_TABLE_POLYHEDRON * ISOSURFACE_TABLE_POLYHEDRON_PTR;

class ISOSURFACE_TABLE {

 protected:
  class ISOSURFACE_TABLE_ENTRY {

  public:
    int num_simplices;
    EDGE_INDEX * simplex_vertex_list;

    ISOSURFACE_TABLE_ENTRY();
    ~ISOSURFACE_TABLE_ENTRY();

    bool Check() const;
  };


 protected:
  ISOSURFACE_TABLE_POLYHEDRON polyhedron;
  ISOSURFACE_TABLE_ENTRY * entry;
  long num_table_entries;

  // max # vertices allowed for table polyhedron
  int max_num_vertices;

  // isosurface table file header
  static const char * const FILE_HEADER;
  static const char * const VERTEX_HEADER;
  static const char * const EDGE_HEADER;
  static const char * const FACET_HEADER;
  static const char * const TABLE_HEADER;

  // old version file header
  static const char * const FILE_HEADER_V1;

  bool is_table_allocated;

  // throw error if table is allocated
  void CheckTableAlloc() const;

 public:
  ISOSURFACE_TABLE(const int d);
  ~ISOSURFACE_TABLE();

  // get functions
  const ISOSURFACE_TABLE_POLYHEDRON & Polyhedron() const
    { return(polyhedron); };
  int Dimension() const { return(polyhedron.Dimension()); };
  int SimplexDimension() const { return(Dimension()-1); };
  int NumVerticesPerSimplex() const { return(SimplexDimension()+1); };
  int NumTableEntries() const { return(num_table_entries); };
  int NumSimplices(const TABLE_INDEX it) const
    { return(entry[it].num_simplices); };
  EDGE_INDEX SimplexVertex
    (const TABLE_INDEX it, const int is, const int iv) const
    // it = table entry index. is = simplex index. iv = vertex index
    // Note: Returns index of edge containing simplex vertex
    { return(entry[it].simplex_vertex_list[is*NumVerticesPerSimplex()+iv]); };
  bool VertexSign(const TABLE_INDEX it, const int iv) const
    { return(polyhedron.VertexSign(it, iv)); };
  int MaxNumVertices() const { return(max_num_vertices); };
  // Note: Even tables for polyhedron of this size are probably impossible 
  //   to compute/store
  const char * const FileHeader() const { return(FILE_HEADER); };
  const char * const VertexHeader() const { return(VERTEX_HEADER); };
  const char * const EdgeHeader() const { return(EDGE_HEADER); };
  const char * const FacetHeader() const { return(FACET_HEADER); };
  const char * const TableHeader() const { return(TABLE_HEADER); };
  const char * const FileHeaderV1() const { return(FILE_HEADER_V1); };
  bool IsTableAllocated() const
    { return(is_table_allocated); };

  // set functions
  void SetNumPolyVertices(const int numv)
    { CheckTableAlloc(); polyhedron.SetNumVertices(numv); };
  void SetNumPolyEdges(const int nume)
    { CheckTableAlloc(); polyhedron.SetNumEdges(nume); };
  void SetNumPolyFacets(const int numf)
    { CheckTableAlloc(); polyhedron.SetNumFacets(numf); };
  void SetPolySize(const int numv, const int nume)
    { SetNumPolyVertices(numv); SetNumPolyEdges(nume); };
  void SetPolySize(const int numv, const int nume, const int numf)
    { SetNumPolyVertices(numv); SetNumPolyEdges(nume); 
      SetNumPolyFacets(numf); };
  void SetPolyVertexCoord(const int iv, const int ic, const int coord)
    { CheckTableAlloc(); polyhedron.SetVertexCoord(iv, ic, coord); };
  // Note: SetNumPolyVertices or SetPolySize must be called before 
  //   SetPolyVertexCoord
  void SetPolyEdge(const int ie, const int iv0, const int iv1)
    { CheckTableAlloc(); polyhedron.SetEdge(ie, iv0, iv1); };
  // Note: SetNumPolyEdges or SetPolySize must be called before SetPolyEdge
  void SetPolyFacetVertex(const int jf, const int iv, const bool in_facet)
    { CheckTableAlloc(); polyhedron.SetFacetVertex(jf, iv, in_facet); };
  // Note: SetPolyNumFacetVertices must be called before SetPolyFacetVertex
  void SetNumSimplices(const TABLE_INDEX it, const int nums);
  void SetSimplexVertex(const TABLE_INDEX it, const int is, 
			const int iv, const EDGE_INDEX ie);
  void Set(const ISOSURFACE_TABLE_POLYHEDRON & new_polyhedron)
    { CheckTableAlloc(); polyhedron = new_polyhedron;
      AllocTable(); };

  // generate polyhedron
  void GenCube(const int cube_dimension)
    { CheckTableAlloc(); polyhedron.GenCube(cube_dimension); AllocTable(); };
  // Note: Cubes of dimension > 4 will have too many vertices
  void GenSimplex(const int simplex_dimension)
    { CheckTableAlloc(); polyhedron.GenSimplex(simplex_dimension); 
      AllocTable(); };

  // allocate memory for isosurface table
  // Must be called after defining polyhedron but before setting simplices
  // Cannot be called twice; Automatically called by GenCube and GenSimplex
  void AllocTable();

  // write isotable
  void WriteText(ostream & out);
  void WriteFileHeaderText(ostream & out);
  void WriteTableText(ostream & out);

  // read isotable
  void ReadText(istream & in);
  void ReadFileHeaderText(istream & in);
  void ReadTableText(istream & in);

  // read/write old versions of isotable file
  void WriteTextV1(ostream & out);
  void ReadTextV1(istream & in);
  void ReadNoHeaderTextV1(istream & in);

  // check functions
  bool CheckDimension(const int d) const;
  bool CheckDimension() const
    { return(CheckDimension(Dimension())); };
  bool CheckTable() const;

};

typedef ISOSURFACE_TABLE * ISOSURFACE_TABLE_PTR;

// routines for generating polyhedra
void generate_prism(const ISOSURFACE_TABLE_POLYHEDRON & base_polyhedron,
		    ISOSURFACE_TABLE_POLYHEDRON & prism);

// error class
class ISOSURFACE_TABLE_ERROR {

 protected:
  char * msg;

 public:
  ISOSURFACE_TABLE_ERROR() { msg = "No error message."; };
  ISOSURFACE_TABLE_ERROR(char * error_msg) 
    { msg = error_msg; };
  char * Msg() const { return(msg); };
};

};

#endif

