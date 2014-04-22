import math;
import matplotlib.pyplot as plt;
import numpy as np;
# ADD-BY-LEETEN 2014/04/19-BEGIN
import scipy.linalg as la;
# ADD-BY-LEETEN 2014/04/19-END
import my_timer;
import incremental_svd as inc_svd;
 
# Parse the 2nd line in the restrt header.
# It should return a tuple of (#atoms, time).
# If the time is not in the file, put 0. 
def parse_restrt_header(header_line):
    header_tokens = header_line.split();
    if( 1 == len(header_tokens) ):
        header_tokens.append("0.0");
    return header_tokens;

# ADD-BY-LEETEN 2014-04-21-BEGIN
def plot_V(fig, V, vmin = -1, vmax = +1):
    plt.figure(fig.number);
    plt.subplot(n_rows, n_cols, subfigi); 
    plt.pcolor(V.T, vmin = vmin, vmax = vmax);
    plt.xlim([0, n_restrts]);
    plt.ylim([0, n_components]);
# ADD-BY-LEETEN 2014-04-21-END
    
restrt_filepath_pattern = "D:/data/mm/Amber_GPU_Benchmark_Suite/PME/JAC_production_NVE/restrt1/restrt1_%d"

step = 100000;
n_restrts = 290;
current_step = 0;  
n_atoms = None;
coord_mat = None;
n_coords_per_atom = 3;
coord_token_width = 12;

def get_subfigure_layout(n_subfigs):
    n_rows = math.floor(math.sqrt(n_subfigs));
    n_cols = math.ceil(n_subfigs/n_rows);
    return int(n_rows), int(n_cols);
    
n_modes = 6;
n_default_dof = 6;
n_components = n_default_dof + n_modes * n_coords_per_atom;
n_training_atoms = 2000;
n_training_coords = n_training_atoms * n_coords_per_atom;

# ADD-BY-LEETEN 2014-04-21-BEGIN
n_update_cols = n_components; 
n_init_cols = n_components;
# ADD-BY-LEETEN 2014-04-21-END

# Compute the pattern of NaN in Amber (*********), and the replacing token (NaN).  
missing_atom_token = "";
for c in range(0, coord_token_width):
    missing_atom_token = missing_atom_token + "*";

new_missing_atom_token = " NaN ";
for c in range(len(new_missing_atom_token), coord_token_width):
    new_missing_atom_token = new_missing_atom_token + " ";

timers = my_timer.create_timers();
my_timer.tic(timers);
for ri in range(0, n_restrts):
    current_step = current_step + step;
    restrt_filepath = restrt_filepath_pattern % current_step;
    with open(restrt_filepath, "rt") as restrt_file:
        # Skip the title.
        restrt_file.readline();
        
        # Read the header.
        restrt_header = restrt_file.readline().strip();
        
        restrt_tokens = parse_restrt_header(restrt_header);
        
        current_n_atoms = int(restrt_tokens[0]);
        
        if( n_atoms is not None and current_n_atoms != n_atoms):
            continue;
        
        n_atoms = current_n_atoms;
        n_coords = n_coords_per_atom * n_atoms;
        
        # np.genfromtxt(restrt_file, delimter=12, missing_values='*');
        coord_buffer = restrt_file.read();
    coord_buffer = "".join(coord_buffer.split("\n"));
    coord_buffer = coord_buffer.replace(missing_atom_token, new_missing_atom_token);
         
    coord_tokens = map(lambda x: coord_buffer[x:x+coord_token_width], range(0, coord_token_width * n_coords, coord_token_width));
        
    coord_list = map(float, coord_tokens);
    coord_vector = np.array([coord_list]).T;
        
    if( coord_mat == None ):
        coord_mat = coord_vector;
    else:
        # MOD-BY-LEETEN 2014-04-21: coord_mat = np.append(coord_mat, coord_vector, axis=1);
        coord_mat = np.hstack([coord_mat, coord_vector]);
        # MOD-BY-LEETEN 2014-04-21-END

coord_mean = np.mean(coord_mat, axis=1);
coord_mean_mat = np.reshape(coord_mean, [len(coord_mean), 1]) * np.ones([1, coord_mat.shape[1]]);
offset_mat = coord_mat - coord_mean_mat;
my_timer.toc(timers);

# TMP: Only take non-WAT atoms.
offset_mat = offset_mat[0:n_training_coords, :];

#############################################################
# incremental SVD.

# # DEL-BY-LEETEN 2014-04-21-BEGIN
# n_update_cols = n_components; 
# n_init_cols = n_components;
# # DEL-BY-LEETEN 2014-04-21-END

my_timer.tic(timers);
# MOD-BY-LEETEN 2014-04-21:    U, s, V = inc_svd.regular_svd_by_pca(offset_mat[:, 0:n_init_cols], k = n_components);
U, s, V = inc_svd.regular_svd(offset_mat[:, 0:n_init_cols], k = n_components);
# MOD-BY-LEETEN 2014-04-21-END
my_timer.toc(timers);

my_timer.tic(timers);
col_base = 0;
col_bases = range(n_init_cols, offset_mat.shape[1], n_update_cols );
n_updates = len(col_bases);
n_rows, n_cols = get_subfigure_layout(n_updates);

# DEL-BY-LEETEN 2014/04/21:    print col_bases;

# MOD-BY-LEETEN 2014/04/19:    plt.figure();
s_fig = plt.figure();
V_fig = plt.figure();
V0_fig = plt.figure();
D_fig = plt.figure();
# MOD-BY-LEETEN 2014/04/19-END
subfigi = 0;
for col_base in col_bases:
    subfigi = subfigi + 1;
    svd_timers = my_timer.create_timers();

    col_end = min(col_base + n_update_cols, offset_mat.shape[1]);
    offset_mat_for_svd = offset_mat[:, 0:col_end]; 

    A = offset_mat[:, col_base:col_end];

    my_timer.tic(svd_timers);
    # MOD-BY-LEETEN 2014/04/19:    U, s, V = inc_svd.addblock_svd_update(U, s, V, A, True);
    U, s, V = inc_svd.addblock(U, s, V, A, True);
    my_timer.toc(svd_timers);

    # Apply regular SVD as the reference.
    my_timer.tic(svd_timers);
    # MOD-BY-LEETEN 2014-04-21:    U0, s0, V0 = inc_svd.regular_svd_by_pca(offset_mat_for_svd, k = n_components);
    U0, s0, V0 = inc_svd.regular_svd(offset_mat_for_svd, k = n_components);
    # MOD-BY-LEETEN 2014-04-21-END
    my_timer.toc(svd_timers);

    my_timer.print_timers("[SVD]", svd_timers);
    
    # Now, compare with the singular values.
    
    # ADD-BY-LEETEN 2014/04/19-BEGIN
    # # MOD-BY-LEETEN 2014-04-21-FROM:
    # plt.figure(V0_fig.number);
    # plt.subplot(n_rows, n_cols, subfigi); 
    # plt.pcolor(V0.T);
    # plt.xlim([0, n_restrts]);
    # plt.ylim([0, n_components]);
    # 
    # V_img = V;
    # for col in range(0, V.shape[1]):
    #     v0 = V0[:, col];
    #     v = V[:, col];
    #     dist = la.norm(v - v0);
    #     negate_dist = la.norm(-v - v0);
    #     if( negate_dist < dist ):
    #         V_img[:, col] = -v;
    #         
    # plt.figure(V_fig.number);
    # plt.subplot(n_rows, n_cols, subfigi); 
    # plt.pcolor(V_img.T);
    # plt.xlim([0, n_restrts]);
    # plt.ylim([0, n_components]);
    # 
    # plt.figure(D_fig.number);
    # plt.subplot(n_rows, n_cols, subfigi); 
    # plt.pcolor(abs((V0 - V_img).T));
    # plt.xlim([0, n_restrts]);
    # plt.ylim([0, n_components]);

    # # MOD-BY-LEETEN 2014-04-21-TO:

    P = V; # np.dot(V, np.diag(s));
    P0 = V0; # np.dot(V0, np.diag(s0));
    for col in range(0, P.shape[1]):
        p0 = P0[:, col];
        p = P[:, col];
        dist = la.norm(p - p0);
        negate_dist = la.norm(-p - p0);
        if( negate_dist < dist ):
            P[:, col] = -p;
    P_diff = np.absolute(P - P0);

    plot_V(V0_fig, P0);            
    plot_V(V_fig, P);            
    plot_V(D_fig, P_diff, vmin = 0, vmax = 1);
    # MOD-BY-LEETEN 2014-04-21-END
    
    plt.figure(s_fig.number);
    # ADD-BY-LEETEN 2014/04/19-END
    plt.subplot(n_rows, n_cols, subfigi); 
    plt.plot(s0);
    plt.plot(s);
    plt.axis('tight');
    # ADD-BY-LEETEN 2014/04/19-BEGIN
    plt.yscale('log');  
    # ADD-BY-LEETEN 2014/04/19-END
    plt.draw();
my_timer.toc(timers);
plt.draw();
plt.show();

my_timer.print_timers("[ALL]", timers);

