/// \file ijkmcube.cxx
/// generate isosurface from scalar field

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

#include <fstream>
#include <iostream>
#include <sstream>

#include "isosurfaceD.h"
#include "ijktable.h"
#include "ijkNrrd.h"
#include "ijkIO.txx"

using namespace ISOSURFACE_D;
using namespace IJKTABLE;

typedef float COORD_TYPE;

typedef enum OUTPUT_FORMAT { OFF, IV };

// global variables
SCALAR_TYPE isovalue = 0.0;
char * filename = NULL;
const char * isotable_directory = 
"/home/4/wenger/programs/ijk/isotable";
OUTPUT_FORMAT output_format = OFF;

// routines
void parse_command_line(int argc, char **argv);
void help(), usage_error();


int main(int argc, char **argv)
{
  try {

  parse_command_line(argc, argv);

  int dimension = 0;
  int * grid_length = NULL;
  SCALAR_TYPE * scalar_grid = NULL;
  Nrrd *nin;

  /* get scalar field data from nrrd file */
  nin = nrrdNew();
  if (nrrdLoad(nin, filename, NULL)) {
    char *err = biffGetDone(NRRD);
    cerr << "Error reading: " << filename << endl;
    cerr << "  Error: " << err << endl;
    exit(35);
  };
  dimension = nin->dim;

  if (dimension < 1) {
    cerr << "Illegal dimension.  Dimension must be at least 1." << endl;
    exit(20);
  };

  grid_length = new int[dimension];
  size_t size[NRRD_DIM_MAX];

  nrrdAxisInfoGet_nva(nin, nrrdAxisInfoSize, size); 

  for (int d = 0; d < dimension; d++)
    grid_length[d] = size[d];

  long num_grid_vert = get_num_grid_vert(dimension, grid_length);

  scalar_grid = new SCALAR_TYPE[num_grid_vert];

  nrrd2scalar(nin, scalar_grid);

  ostringstream isotable_pathname;
  ostringstream isotable_filename;

  isotable_filename.str("");
  isotable_pathname.str("");
  isotable_filename << "isosurface_table.c." << dimension << ends;
  isotable_pathname << isotable_directory << "/" 
		    << isotable_filename.str() << ends;

  ifstream isotable_file(isotable_pathname.str().c_str(), ios::in);
  if (!isotable_file) {
    isotable_file.clear();
    isotable_file.open(isotable_filename.str().c_str(), ios::in);
  };

  if (!isotable_file) {
    cerr << "Unable to obtain isosurface table file "
	 << isotable_pathname.str() << "." << endl;
    exit(30);
  };

  ISOSURFACE_TABLE isotable(dimension);

  try {
    isotable.ReadText(isotable_file);

    if (isotable.Dimension() != dimension) {
      cerr << "Error.  Isotable dimension does not equal volume dimension."
	   << endl;
      cerr << "  Isotable dimension = " << isotable.Dimension()
	   << ".  Volume dimension = " << dimension << "." << endl;
      exit(40);
    };
  }
  catch (ISOSURFACE_TABLE_ERROR error) {
    cerr << "Error reading isosurface table from " 
	 << isotable_pathname.str() << "." << endl;
    cerr << error.Msg() << endl;
    cerr << "Exiting." << endl;
    exit(35);
  };

  int num_intersected_edges = 0;
  int num_simplices = 0;
  int * intersected_edge_endpoint = NULL;
  int * simplex_vert = NULL;

  isosurfaceD(dimension, isotable, scalar_grid, grid_length, isovalue, 
	      intersected_edge_endpoint, num_intersected_edges,
	      simplex_vert, num_simplices);


  int numv = num_intersected_edges;
  COORD_TYPE * coord = new COORD_TYPE[dimension*numv];
  int * coord0 = new int[dimension];
  int * coord1 = new int[dimension];

  for (int iv = 0; iv < numv; iv++) {
    int v0 = intersected_edge_endpoint[2*iv];
    int v1 = intersected_edge_endpoint[2*iv+1];
    SCALAR_TYPE s0 = scalar_grid[v0];
    SCALAR_TYPE s1 = scalar_grid[v1];
    COORD_TYPE w0, w1;
    double s_diff = s1 - s0;
    const double EPSILON = 0.00001;
    if (s_diff > EPSILON || s_diff < -EPSILON) { 
      w0 = (s1 - isovalue) / s_diff;
      w1 = (isovalue - s0) / s_diff;
    }
    else {
      // arbitrarily set weights to 0.5
      w0 = w1 = 0.5;
    };

    get_coord(v0, dimension, grid_length, coord0);
    get_coord(v1, dimension, grid_length, coord1);

    for (int d = 0; d < dimension; d++)
      coord[iv*dimension+d] = w0*coord0[d] + w1*coord1[d];
  }
  delete [] coord0;
  delete [] coord1;

  try {

    IJKIO::ERROR error_mcube("ijkmcube");

    switch (output_format) {

    case OFF:
      IJKIO::ijkoutOFF(dimension, coord, numv, simplex_vert, num_simplices);
      break;

    case IV:
      if (dimension == 3) {
	IJKIO::ijkoutIV(dimension, coord, numv, simplex_vert, num_simplices);
      }
      else throw error_mcube("Illegal dimension. OpenInventor format is only for dimension 3.");
      break;

    default:
      throw error_mcube("Illegal output format.");
      break;
    }
  }
  catch(IJKIO::ERROR error) {
    cerr << "Error in procedure: " << error.procedure << endl;
    cerr << "  " << error.msg1 << endl;
  }

  delete [] intersected_edge_endpoint;
  delete [] simplex_vert;
  delete [] coord;
  delete [] scalar_grid;
  delete [] grid_length;

  nrrdNuke(nin);

  } 
  catch (ISOSURFACE_TABLE_ERROR error) {
    cerr << "Programming error." << endl;
    cerr << error.Msg() << endl;
    cerr << "Exiting." << endl;
    exit(60);
  };

}

void parse_command_line(int argc, char **argv)
{
  int iarg = 1;
  while (iarg < argc && argv[iarg][0] == '-') {
    if (strcmp(argv[iarg], "-dir") == 0) {
      iarg++;
      if (iarg >= argc) usage_error();

      isotable_directory = argv[iarg];
    }
    else if (strcmp(argv[iarg], "-off") == 0) {
      output_format = OFF;
    }
    else if (strcmp(argv[iarg], "-iv") == 0) {
      output_format = IV;
    }
    else if (strcmp(argv[iarg], "-help") == 0) {
      help();
    }
    else
      usage_error();

    iarg++;
  };

  if ((argc - iarg) != 2) usage_error();

  sscanf(argv[iarg], "%f", &isovalue);
  filename = argv[iarg+1];
}

void usage_msg()
{
  cerr << "Usage: ijkmcube [-help] [-off|-iv] [-dir {isotable_directory}] {isovalue} {filename}" 
       << endl;
}

void usage_error()
{
  usage_msg();
  exit(10);
}

void help()
{
  usage_msg();
  cerr << endl;
  cerr << "ijkmcube - Marching cubes isosurface generation algorithm." << endl;
  cerr << endl;
  cerr << "Options:" << endl;

  cerr << "  -off: Output in geomview OFF format. (Default.)" << endl;
  cerr << "  -iv: Output in OpenInventor .iv format." << endl;
  cerr << "  -dir {isotable_directory}: Directory containing appropriate isosurface table." << endl;
  cerr << "  -help: Print this help message." << endl;
  exit(20);
}
