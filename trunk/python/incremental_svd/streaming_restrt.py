import math;
import matplotlib.pyplot as plt;
import numpy as np;
import scipy.linalg as la;
import sys;
import traceback;
import my_md;
import my_timer;
import incremental_svd as inc_svd;
from optparse import OptionParser;

def get_subfigure_layout(n_subfigs):
    n_rows = math.floor(math.sqrt(n_subfigs));
    n_cols = math.ceil(n_subfigs/n_rows);
    return int(n_rows), int(n_cols);

restrt_filepath_pattern = "D:/data/mm/Amber_GPU_Benchmark_Suite/PME/JAC_production_NVE/restrt1/restrt1_%d"

# # TEST-MOD-FROML
# step = 100000;
# n_restrts = 290;
# # TO:
step = 20000;
step_begin = 0;
n_restrts = 1000;
# # TEST-MOD-END

# # TEST-MOD-FROM:
# n_modes = 6;
# n_default_dof = 6;
# n_components = n_default_dof + n_modes * n_coords_per_atom;
# # TO:
n_components = 200;
with_ground_truth = False;
# TEST-MOD-END
n_training_atoms = 2000;

n_update_cols = 200; 
n_init_cols = 200;

parser = OptionParser()
parser.add_option("--restrtFilepathPrefix", default=restrt_filepath_pattern);
parser.add_option("--nrOfRestrts",          default=n_restrts, type=int);
parser.add_option("--restrtBegin",            default=step_begin, type=int);
parser.add_option("--restrtStep",         default=step, type=int);
parser.add_option("--withGroundTruth",      action="store_true", dest = "withGroundTruth");
parser.add_option("--nrOfComponents",       default=n_components, type=int);
parser.add_option("--nrOfUpdateColumns",    default=n_update_cols, type=int);
parser.add_option("--nrOfInitColumns",      default=n_init_cols, type=int);
parser.add_option("--nrOfTrainingAtoms",    default=n_training_atoms, type=int);
(options, args) = parser.parse_args()

restrt_filepath_pattern = options.restrtFilepathPrefix;
n_restrts = options.nrOfRestrts;
step_begin = options.restrtBegin;
step = options.restrtStep;
with_ground_truth = options.withGroundTruth;
n_components = options.nrOfComponents;
n_update_cols = options.nrOfUpdateColumns;
n_init_cols = options.nrOfInitColumns;
n_training_atoms = options.nrOfTrainingAtoms;

n_atoms = None;
coord_mat = None;
n_coords_per_atom = my_md.get_n_coords_per_atom();
n_training_coords = n_training_atoms * n_coords_per_atom;

col_bases = range(n_init_cols, n_restrts, n_update_cols );
n_updates = len(col_bases);
n_rows, n_cols = get_subfigure_layout(n_updates);

s_fig = plt.figure();
subfigi = 0;

# Frobenius norms.
Fs = [];
current_step = step_begin;
for ri in range(0, n_restrts):
    current_step = current_step + step;
    restrt_filepath = restrt_filepath_pattern % current_step;
    coord_vector = my_md.read_restrt_coord(restrt_filepath);
    
    coord_vector = coord_vector[0:n_training_coords];
    if( coord_mat == None ):
        coord_mat = coord_vector;
    else:
        coord_mat = np.hstack([coord_mat, coord_vector]);
        
    if( coord_mat.shape[1] == n_init_cols ):
        U, s, V = inc_svd.regular_svd(coord_mat, k = n_components);
        coord_sum = np.sum(coord_mat, axis=1);
        
        A = coord_mat; 
        coord_mean = np.mean(A, axis=1);  
        U_offset, s_offset, V_offset = inc_svd.update(U, s, V, -coord_mean.reshape(-1, 1), np.ones([V.shape[0], 1]), True);
        
        # Check how well the new U_offset can represent the offset of A
        A_offset = A - np.dot(coord_mean.reshape(-1, 1), np.ones([1, A.shape[1]]));
        A_diff_from_U = A_offset - np.dot(U, np.dot(U.T, A_offset));
        f = la.norm(A_diff_from_U);
        Fs.append(np.sqrt(f * f /np.prod(A.shape)));
        

    if( coord_mat.shape[1] > n_init_cols ):
        A = None;
        if( 0 == (coord_mat.shape[1] - n_init_cols) % n_update_cols ):
            A = coord_mat[:, -n_update_cols::];
        elif ( coord_mat.shape[1] == n_restrts ):
            A = coord_mat[:, n_init_cols + n_update_cols * int((coord_mat.shape[1] - n_init_cols)/n_update_cols)::];

        if( A is not None ):
            timers = my_timer.create_timers();
            my_timer.tic(timers);
            U, s, V = inc_svd.addblock(U, s, V, A, True);
            my_timer.toc(timers);
            
            my_timer.tic(timers);
            coord_sum = coord_sum + np.sum(A, axis=1);
            coord_mean = coord_sum / float(coord_mat.shape[1]);  
            U_offset, s_offset, V_offset = inc_svd.update(U, s, V, -coord_mean.reshape(-1, 1), np.ones([V.shape[0], 1]), True);
            my_timer.toc(timers);
            
            # Check how well the new U_offset can represent the offset of A
            my_timer.tic(timers);
            A_offset = A - np.dot(coord_mean.reshape(-1, 1), np.ones([1, A.shape[1]]));
            A_diff_from_U = A_offset - np.dot(U, np.dot(U.T, A_offset));
            f = la.norm(A_diff_from_U);
            Fs.append(np.sqrt(f * f /np.prod(A.shape)));
            my_timer.toc(timers);
            
            subfigi = subfigi + 1; 
            plt.figure(s_fig.number);
            plt.subplot(n_rows, n_cols, subfigi); 
            plt.plot(s_offset);
            
            if( with_ground_truth ):
                my_timer.tic(timers);
                try:
                    offset_mat = coord_mat - np.dot(coord_mean.reshape(-1, 1), np.ones([1, coord_mat.shape[1]])); 
                    U0, s0, V0 = inc_svd.regular_svd(offset_mat, k=n_components);
                    plt.plot(s0);
                except Exception, e:
                    traceback.print_exc(file=sys.stdout);
                my_timer.toc(timers);
                
            plt.axis('tight');
            plt.yscale('log');  
            
            my_timer.print_timers("[SVD]", timers);
            
plt.figure();
plt.plot(Fs);
plt.axis('tight');

plt.show();