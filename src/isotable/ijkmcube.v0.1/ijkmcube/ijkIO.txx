// ijk output templates
// output OpenInventor .iv file
// output Geomview file

#ifndef _IJKIO_
#define _IJKIO_

#include <iostream>
#include <string>

namespace IJKIO {

  //******************************************
  // Error handler
  //******************************************

  class ERROR {

  public:
    std::string procedure;       // procedure name
    std::string msg1;            // error message 1 (main error message)
    std::string msg2;            // error message 2 (hints on cause)
    int num;                     // error number

  protected:
    void Init()
    { procedure = msg1 = msg2 = ""; num = 0; }

  public:
    ERROR(const std::string proc)
    { Init(); this->procedure = proc; };

    ERROR & operator ()(const std::string msg1)
    { this->msg1 = msg1; return(*this); };

    void operator ()(const std::string msg1, const std::string msg2)
    { this->msg1 = msg1; this->msg2 = msg2; };

  };

  //******************************************
  // Output OpenInventor .iv file
  //******************************************

  template <class T> void ijkoutIV
  (std::ostream & out, const int dim, const T * coord, const int numv,
   const int * tri, const int numt)
    // output OpenInventor .iv file
    // out = output stream
    // dim = dimension.  Must be 3.
    // coord[3*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // tri[3*j+k] = k'th vertex index of triangle j (k < dim)
    // numt = number of triangles
  {
    ERROR error("ijkoutIV");

    if (dim != 3)
      throw error("Illegal dimension.  OpenInventor files are only for dimension 3.");

    out << "#Inventor V2.1 ascii" << endl;
    out << endl;

    out << "Separator {" << endl;

    // set vertex ordering to clockwise to turn on two-sided lighting
    out << "  ShapeHints {" << endl;
    out << "    vertexOrdering CLOCKWISE" << endl;
    out << "  }" << endl;
    out << endl;

    out << "  IndexedFaceSet {" << endl;

    out << "    vertexProperty VertexProperty {" << endl;
    out << "      vertex [" << endl;

    // output vertex coordinates
    out << endl << "# vertex coordinates" << endl;
    for (int i = 0; i < numv; i++) {
      for (int d = 0; d < dim; d++) {
	out << coord[dim*i+d];
	if (d < dim-1) { out << " "; }
	else {	
	  if (i < numv-1) { out << "," << endl; };
	};
      };
    };

    out << " ]" << endl;
    out << "    }" << endl;

    out << "    coordIndex [" << endl;
    // output triangle verices
    out << endl << "# triangle vertices" << endl;
    for (int it = 0; it < numt; it++) {
      for (int d = 0; d < dim; d++) {
	out << tri[dim*it+d] << ",";
      };
      out << "-1";
      if (it < numt-1) {
	out << "," << endl;
      };
    };
    out << " ]" << endl;

    out << "  }" << endl;

    out << "}" << endl;
  }

  template <class T> void ijkoutIV
  (const int dim, const T * coord, const int numv,
   const int * tri, const int numt)
    // output OpenInventor .iv format to standard outpu
    // dim = dimension.  Must be 3.
    // coord[3*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // tri[3*j+k] = k'th vertex index of triangle j (k < dim)
    // numt = number of triangles
  {
    ijkoutIV(cout, dim, coord, numv, tri, numt); 
  }

  //******************************************
  // Geomview OFF file
  //******************************************

  template <class T> void ijkoutOFF
  (std::ostream & out, const int dim, const T * coord, const int numv,
   const int * simplex_vert, const int nums)
    // output Geomview OFF files
    // out = output stream
    // dim = dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
  {
    ERROR error("ijkoutOFF");

    if (dim == 3) { out << "OFF" << endl; }
    else if (dim == 4) { out << "4OFF" << endl;}
    else {
      out << "nOFF" << endl;
      out << dim << endl;
    };

    out << numv << " " << nums << " " << 0 << endl;

    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) {
	out << coord[iv*dim + d];
	if (d < dim-1) { out << " "; }
	else { out << endl; };
      }
    };
    out << endl;

    for (int is = 0; is < nums; is++) {
      out << dim << " ";
      for (int d = 0; d < dim; d++) {
	out << simplex_vert[is*dim + d];
	if (d < dim-1) { out << " "; }
	else { out << endl; };
      };
    };

  }

  template <class T> void ijkoutOFF
  (const int dim, const T * coord, const int numv,
   const int * simplex_vert, const int nums)
    // output Geomview OFF format to standard output
    // dim = dimension.
    // coord[3*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[3*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
  {
    ijkoutOFF(cout, dim, coord, numv, simplex_vert, nums);
  }

  template <class T> void ijkinOFF
  (std::istream & in, int & dim, T * & coord, int & numv,
   int * & simplex_vert, int & nums)
    // input Geomview OFF files
    // in = input stream
    // dim = dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
  {
    std::string header_line;
    int nume;

    ERROR error("ijkinOFF");

    coord = NULL;
    simplex_vert = NULL;

    // read input
    header_line = "";
    while (header_line == "" && !in.eof())
      in >> header_line;

    if (header_line == "OFF") {
      dim = 3;
    }
    else if (header_line == "4OFF") {
      dim = 4;
    }
    else if (header_line == "nOFF") {
      in >> dim;
    }
    else {
      std::string errmsg =  
	"Illegal Geomview .off file header: " + header_line;
      throw error(errmsg);
    };

    if (dim < 1)
      throw error("Dimension must be at least 1.");

    in >> numv;
    in >> nums;
    in >> nume;

    coord = new T[numv*dim];
    for (int iv = 0; iv < numv; iv++) {
      for (int d = 0; d < dim; d++) 
	in >> vertex_coord[iv*dim + d];
    };

    simplex_vert = new int[nums*dim];
    for (int is = 0; is < nums; is++) {

      int num_simplex_vert = 0;
      in >> num_simplex_vert;
      if (num_simplex_vert != dim) {
	delete [] coord;
	delete [] simplex_vert;
	coord = NULL;
	simplex_vert = NULL;
	throw error("Wrong number of vertices in geomview polytope list.");
      }
      for (int d = 0; d < dim; d++)
	in >> simplex_vert[is*dim + d];
    };
  }

  template <class T> void ijkinOFF
  (int & dim, T * & coord, int & numv, int * & simplex_vert, int & nums)
    // input Geomview OFF format to standard input
    // in = input stream
    // dim = dimension
    // coord[dim*i+k] = k'th coordinate of vertex j  (k < dim)
    // numv = number of vertices 
    // simplex_vert[dim*j+k] = k'th vertex index of simplex j
    // nums = number of simplices
  {
    ijkinOFF(cin, dim, coord, numv, simplex_vert, nums);
  }

}

#endif
