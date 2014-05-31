# This is the python porting of D. Wingate's MATLAB script of incremental SVD. 
# The comments and logic flow follow D. Wingate's implementation.
#
# The algorithm is based on Matthew Brand's incremetnal SVD algorithm.
# "Fast low-rank modifications of the thin singular value decomposition".
#
# Teng-Yok Lee, 2014/04/19.

import sys;
import numpy as np;
import scipy.linalg as la;
import scipy.sparse as sparse;
import scipy.sparse.linalg as sparsela; 

# If the matrix is very thin, using PCA will be fast.
def regular_svd_by_pca(K, k=0):
    
    K_size = K.shape;
    if( K_size[0] < K_size[1] ):
        K_squared = np.dot(K, K.T);        
        tsp, tUp = la.eig(K_squared);
    else:
        K_squared = np.dot(K.T, K);
        tsp, tVp = la.eig(K_squared);
    # As la.eig returns complex number, use its absolute value.
    tsp = abs(tsp);
    tsp = np.sqrt(tsp);
    n_pos_sigs = sum(tsp > 0);
    tSp = np.diag(map(lambda s: 1.0/s, tsp[0:n_pos_sigs]));
    if( K_size[0] < K_size[1] ):
        tVp = np.dot(K.T, tUp);
        tVp[:, 0:n_pos_sigs] = np.dot(tVp[:, 0:n_pos_sigs], tSp);
    else:
        tUp = np.dot(K, tVp);
        tUp[:, 0:n_pos_sigs] = np.dot(tUp[:, 0:n_pos_sigs], tSp);
    if( 0 < k and k < min(K_size) ):
        tUp = tUp[:, 0:k];
        tVp = tVp[:, 0:k];
        tsp = tsp[0:k];
    return tUp, tsp, tVp;

def orth(A, threshold=1e-6):
    q, r = np.linalg.qr(A, mode='reduced');
    return q;

# A wrapper of svd.
# Otherwise, scipy.sparse.linals's svd cannot handle the case that the rank is equal to the min(#rows, #col).  
def regular_svd(A, k=0):
    if( k <= 0 or k >= min(A.shape) ):
        U, s, VT = la.svd(A, full_matrices=False);
    else:    
        U, s, VT = sparsela.svds(sparse.bsr_matrix(A), k = k);
        U = U[:,::-1];
        s = s[::-1];
        VT = VT[::-1, :];
        
    V = VT.T;
    
    # If the actual rank is smaller than k, part of the metrics will be nan.
    if( not np.all(np.isfinite(s)) ):
        print >> sys.stderr, "[WARN] regular_svd: s has non-finite numbers."
        finite_indices = np.isfinite(s);
        U = U[:, finite_indices]; 
        s = s[finite_indices]; 
        V = V[:, finite_indices];
    return U, s, V; 
    

#! Given X = USV', compute SVD for X + A * B'.
#
# Based on http://web.mit.edu/~wingated/www/scripts/svd_update.m.
#
def update( U, s, V, A, B, force_orth ):
    current_rank = U.shape[1];
    
    # P: An orthogonal basis of the column-space of (I-UU')A, 
    # the component of "A" that is orthogonal to U.
    m = np.dot(U.T, A);
    p = A - np.dot(U, m);

    P = orth( p );
    
    # Pad with zeros if p does not have full rank.
    P = np.pad(P, ((0,0), (0, p.shape[1] - P.shape[1])), 'constant', constant_values=(0, 0));
    
    Ra = np.dot(P.T, p);

    #: Q: An orthogonal basis of the column-space of (I-VV')b.
    n = np.dot(V.T, B);
    q = B - np.dot(V, n);

    Q = orth(q);
    
    # Pad with zeros if q does not have full rank.
    Q = np.pad(Q, ((0,0), (0, q.shape[1] - Q.shape[1])), 'constant', constant_values=(0, 0));

    Rb = np.dot(Q.T, q);

    K = np.pad(np.diag(s), ((0, m.shape[1]), (0, m.shape[1])), 'constant', constant_values=(0, 0)) + np.dot(np.vstack([m, Ra]), np.hstack([n.T, Rb.T]));

    tUp, tsp, tVp = regular_svd(K, k = current_rank);
    
    # update our matrices.

    sp = tsp;
    Up = np.dot(np.hstack([U, P]), tUp);
    Vp = np.dot(np.hstack([V, Q]), tVp);

    # The above rotations may not preserve orthogonality, so we explicitly
    # deal with that via a QR plus another SVD.  In a long loop, you may
    # want to force orthogonality every so often.

    if ( force_orth ):
        Up, sp, Vp = re_orth(Up, sp, Vp);
  
    return [Up, sp, Vp];

#! Given X = USV', compute SVD when appending new columns A to X
#
# Based on http://web.mit.edu/~wingated/www/scripts/addblock_svd_update.m
#

def addblock( U, s, V, A, force_orth ):
    
    current_rank = U.shape[1];

    # P is an orthogonal basis of the column-space
    # of (I-UU')a, which is the component of "a" that is
    # orthogonal to U.
    m = np.dot(U.T, A);
    p = A - np.dot(U, m);

    P = orth( p );
    if( P.size > 0 ):
        P = np.pad(P, ((0,0), (0, p.shape[1] - P.shape[1])), 'constant', constant_values=(0, 0));
        Ra = np.dot(P.T, p);

        #  
        # Diagonalize K, maintaining rank
        #
    
        z = np.zeros( m.shape );
        K = np.vstack([np.hstack([np.diag(s), m]), np.hstack([z.T, Ra])]);
    else:
        K = np.hstack([np.diag(s), m]);
    
    tUp, tsp, tVp = regular_svd(K, k = current_rank);

    #
    # Now update our matrices!
    #

    sp = tsp;

    # this may not preserve orthogonality over
    # many repetitions.  See below.
    Up = np.dot(np.hstack([U, P]), tUp);   

    # Exploit structure to compute this fast: Vp = [ V Q ] * tVp;
    if( 0 != V.size ):
        Vp = np.dot(V, tVp[0:current_rank, : ]);
    else:
        Vp = np.zeros([0, 0]);
    Vp = np.vstack([
                    Vp, 
                    tVp[current_rank:tVp.shape[0], : ]
                    ]);

    # The above rotations may not preserve orthogonality, so we explicitly
    # deal with that via a QR plus another SVD.  In a long loop, you may
    # want to force orthogonality every so often.

    if ( force_orth ):
        Up, sp, Vp = re_orth(Up, sp, Vp);

    return [Up, sp, Vp];

# TODO: Port the rank one update too.
# http://web.mit.edu/~wingated/www/scripts/rank_one_svd_update.m

def my_addblock( U, s, V, A, force_orth ):
    new_n_cols = A.shape[1];
    V_padded = np.pad(V, ((0, new_n_cols), (0, 0)), 'constant', constant_values=(0, 0));
    B = np.pad(np.eye(new_n_cols, new_n_cols), ((V.shape[0], 0), (0, 0)), 'constant', constant_values=(0, 0));
    return update(U, s, V_padded, A, B, force_orth);

#! Given X = USV', compute SVD for X + a * b'.
#
# Based on http://web.mit.edu/~wingated/www/scripts/rank_one_svd_update.m
#
def rank_one_update( U, s, V, a, b, force_orth ):

    current_rank = U.shape[1];

    K = np.diag(s);
    
    m = np.dot(U.T, a);
    p = a - np.dot(U, m);
    # P = p / la.norm(p);
    # Ra = np.dot(P.T, p);
    p_norm = la.norm(p); 
    if (p_norm  > 1e-6 ):
        P = p / p_norm;
        Ra = np.dot(P.T, p);
        K = np.pad(K, ((0, 1), (0, 0)), 'constant', constant_values=(0,0));
        m = np.append(m, Ra);
    else:
        P = np.zeros([0, 0]);
        
    # % Q is an orthogonal basis of the column-space
    # % of (I-VV')b.
    n = np.dot(V.T, b);
    q = b - np.dot(V, n);
    # Q = q / la.norm(q);
    # Rb = np.dot(Q.T, q);
    q_norm = la.norm(q); 
    if (q_norm  > 1e-6 ):
        Q = q / q_norm;
        Rb = np.dot(Q.T, q);
        K = np.pad(K, ((0, 0), (0, 1)), 'constant', constant_values=(0,0));
        n = np.append(n, Rb);
    else:
        Q = np.zeros([0, 0]);

    #  %
    #  % Diagonalize K, maintaining rank
    #  %
    K = K + np.dot(m.reshape(-1, 1), n.reshape(-1, 1).T);
    tUp, tsp, tVp = regular_svd( K, current_rank );

    # %
    # % Now update our matrices!
    # %
  
    sp = tsp;
    if( P.size > 0 ):
        Up = np.dot(np.hstack([U, P.reshape(-1, 1)]), tUp);
    else:
        Up = np.dot(U, tUp);
        
    if( Q.size > 0 ):
        Vp = np.dot(np.hstack([V, Q.reshape(-1, 1)]), tVp);
    else:
        Vp = np.dot(V, tVp);
    
    # % The above rotations may not preserve orthogonality, so we explicitly
    # % deal with that via a QR plus another SVD.  In a long loop, you may
    # % want to force orthogonality every so often.

    if ( force_orth ):
        Up, sp, Vp = re_orth(Up, sp, Vp);
        
    return [Up, sp, Vp];

def re_orth(Up, sp, Vp):
    UQ, UR = np.linalg.qr(Up, mode='reduced');
    VQ, VR = np.linalg.qr(Vp, mode='reduced');
    
    [tUp, tsp, tVp] = regular_svd( np.dot(np.dot(UR, np.diag(sp)), VR.T) );
    return np.dot(UQ, tUp), tsp, np.dot(VQ, tVp);
    
def addvector( U, s, V, a, force_orth ):
    V_padded = np.pad(V, ((0, 1), (0, 0)), 'constant', constant_values=(0, 0));
    b = np.pad(np.eye(1, 1), ((V.shape[0], 0), (0, 0)), 'constant', constant_values=(0, 0));
    return rank_one_update(U, s, V_padded, a, b, force_orth);    

#! Apply hotelling transform to cancel translation/rotation of the entire system.
#
def apply_hotelling_transform(coord_mat):
    coord_mean = np.mean(coord_mat, axis=1);
    coord_mean_mat = np.reshape(coord_mean, [len(coord_mean), 1]) * np.ones([1, coord_mat.shape[1]]);
    offset_mat = coord_mat - coord_mean_mat;
    U, s, V = regular_svd(offset_mat);
    return np.dot(U.T, offset_mat);
    
