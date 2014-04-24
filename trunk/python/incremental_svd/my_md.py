import os;
import numpy as np;

n_coords_per_atom = 3;
coord_token_width = 12;
  
# Compute the pattern of NaN in Amber (*********), and the replacing token (NaN).
missing_atom_token = "";
for c in range(0, coord_token_width):
    missing_atom_token = missing_atom_token + "*";

new_missing_atom_token = " NaN ";
for c in range(len(new_missing_atom_token), coord_token_width):
    new_missing_atom_token = new_missing_atom_token + " ";

def parse_restrt_header(header_line):
    header_tokens = header_line.split();
    if( 1 == len(header_tokens) ):
        header_tokens.append("0.0");
    return header_tokens;

def get_n_coords_per_atom():
    return n_coords_per_atom;

def read_restrt_coord(restrt_filepath):

    if( os.path.exists(restrt_filepath) == False ):
        return None;
    
    with open(restrt_filepath, "rt") as restrt_file:
        # Skip the title.
        restrt_file.readline();
        
        # Read the header.
        restrt_header = restrt_file.readline().strip();
        
        restrt_tokens = parse_restrt_header(restrt_header);
        
        n_atoms = int(restrt_tokens[0]);
        n_coords = n_coords_per_atom * n_atoms;
        
        # np.genfromtxt(restrt_file, delimter=12, missing_values='*');
        coord_buffer = restrt_file.read();
    coord_buffer = "".join(coord_buffer.split("\n"));
    coord_buffer = coord_buffer.replace(missing_atom_token, new_missing_atom_token);
         
    coord_tokens = map(lambda x: coord_buffer[x:x+coord_token_width], range(0, coord_token_width * n_coords, coord_token_width));
        
    coord_list = map(float, coord_tokens);
    coord_vector = np.array([coord_list]).T;

    return coord_vector;