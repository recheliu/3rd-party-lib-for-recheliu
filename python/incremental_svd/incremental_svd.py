# This is the python porting of D. Wingate's MATLAB script of incremental SVD. 
# The comments and logic flow follow D. Wingate's implementation.
#
# The algorithm is based on Matthew Brand's incremetnal SVD algorithm.
# "Fast low-rank modifications of the thin singular value decomposition".
#
# Teng-Yok Lee, 2014/04/19.

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

# ADD-BY-LEETEN 2014/04/19-BEGIN
# # MOD-BY-LEETEN 2014/04/20-FROM:
# def orth(A):
#     U, s, V = regular_svd_by_pca(A, k=min(A.shape));
#     return U;
# # MOD-BY-LEETEN 2014/04/20-TO:
def orth(A, threshold=1e-6):
    U, s, V = regular_svd(A, k=min(A.shape));
    return U[:, abs(s) > threshold];
# # MOD-BY-LEETEN 2014/04/20-END
# ADD-BY-LEETEN 2014/04/19-END

# A wrapper of svd.
# Otherwise, scipy.sparse.linals's svd cannot handle the case that the rank is equal to the min(#rows, #col).  
def regular_svd(A, k=0):
    # # MODBY-LEETEN 2014/04/20-FROM:
    # if( k <= 0 or k >= min(A.shape) ):
    #     return la.svd(A, full_matrices=False);
    # else:    
    #     return sparsela.svds(sparse.bsr_matrix(A), k = k);
    # # MODBY-LEETEN 2014/04/20-TO:
    if( k <= 0 or k >= min(A.shape) ):
        U, s, VT = la.svd(A, full_matrices=False);
    else:    
        U, s, VT = sparsela.svds(sparse.bsr_matrix(A), k = k);
        U = U[:,::-1];
        s = s[::-1];
        VT = VT[::-1, :];
    return U, s, VT.T;
    # # MODBY-LEETEN 2014/04/20-END
    

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

    # MOD-BY-LEETEN 2014/04/19:    P = la.orth( p );
    P = orth( p );
    # MOD-BY-LEETEN 2014/04/19-END
    # Pad with zeros if p does not have full rank.
    # MOD-BY-LEETEN 2014/04/19:    P = np.hstack([P, np.zeros(P.shape[0], p.shape[1] - P.shape[1])]);
    P = np.pad(P, ((0,0), (0, p.shape[1] - P.shape[1])), 'constant', constant_values=(0, 0));
    # MOD-BY-LEETEN 2014/04/19-END

    Ra = np.dot(P.T, p);

    #: Q: An orthogonal basis of the column-space of (I-VV')b.
    n = np.dot(V.T, B);
    q = B - np.dot(V, n);

    # MOD-BY-LEETEN 2014/04/19:    Q = la.orth( q );
    Q = orth(q);
    # MOD-BY-LEETEN 2014/04/19-END
    # Pad with zeros if q does not have full rank.
    # MOD-BY-LEETEN 2014/04/19:    Q = np.hstack([Q, np.zeros(Q.shape[0], q.shape[1] - Q.shape[1])]);
    # MOD-BY-LEETEN 2014/04/21:    Q = np.pad(Q, ((0,0), (0, q.shape[1] - q.shape[1])), 'constant', constant_values=(0, 0));
    Q = np.pad(Q, ((0,0), (0, q.shape[1] - Q.shape[1])), 'constant', constant_values=(0, 0));
    # MOD-BY-LEETEN 2014/04/21-END
    # MOD-BY-LEETEN 2014/04/19-END 

    Rb = np.dot(Q.T, q);

    # # MOD-BY-LEETEN 2014/04/21-FROM: 
    # # Diagonalize K, maintaining rank.
    # z = np.zeros( m.shape );
    # z2 = np.zeros( [m.shape[1], m.shape[1]] );
    # 
    # # MOD-BY-LEETEN 2014/04/19:    K = np.vstack([np.hstack([np.diag(s), z]), np.hstack([z.T, z2])]) + np.vstack([m, Ra]) * np.hstack(n, Rb);
    # K = np.vstack([np.hstack([np.diag(s), z]), np.hstack([z.T, z2])]) + np.dot(np.vstack([m, Ra]), np.hstack([n.T, Rb.T]));
    # # MOD-BY-LEETEN 2014/04/19-END
    # # MOD-BY-LEETEN 2014/04/21-TO:
    K = np.pad(np.diag(s), ((0, m.shape[1]), (0, m.shape[1])), 'constant', constant_values=(0, 0)) + np.dot(np.vstack([m, Ra]), np.hstack([n.T, Rb.T]));
    # # MOD-BY-LEETEN 2014/04/21-END

    # tUp, tsp, tVpT = sparsela.svds(sparse.bsr_matrix(K), k=current_rank );
    # tUp = tUp[:, ::-1];
    # tSp = np.diag(tsp[::-1]); 
    # tVp = tVpT[::-1, :].T;
    # MOD-BY-LEETEN 2014/04/20:    tUp, tsp, tVp = regular_svd_by_pca(K, k = current_rank);
    tUp, tsp, tVp = regular_svd(K, k = current_rank);
    # MOD-BY-LEETEN 2014/04/20-END
  
    # update our matrices.

    sp = tsp;
    Up = np.dot(np.hstack([U, P]), tUp);
    Vp = np.dot(np.hstack([V, Q]), tVp);

    # The above rotations may not preserve orthogonality, so we explicitly
    # deal with that via a QR plus another SVD.  In a long loop, you may
    # want to force orthogonality every so often.

    if ( force_orth ):
        # # MOD-BY-LEETEN 2014/04/19-FROM:
        # UQ, UR = la.qr( Up, model='economic' );
        # VQ, VR = la.qr( Vp, model='economic' );
        # # TO:
        # # MOD-BY-LEETEN 2014/04/21-FROM:
        # UQ, UR = la.qr( Up, mode='economic' );
        # VQ, VR = la.qr( Vp, mode='economic' );
        # # # MOD-BY-LEETEN 2014/04/19-END
        # [tUp, tsp, tVp] = la.svd( np.dot(np.dot(UR, np.diag(sp)), VR.T) );
        # Up = np.dot(UQ, tUp);
        # Vp = np.dot(VQ, tVp);
        # sp = tsp;
        # # MOD-BY-LEETEN 2014/04/21-TO:
        Up, sp, Vp = re_orth(Up, sp, Vp);
        # # MOD-BY-LEETEN 2014/04/21-END
  
    return [Up, sp, Vp];

#! Given X = USV', compute SVD when appending new columns A to X
#
# Based on http://web.mit.edu/~wingated/www/scripts/addblock_svd_update.m
#

# MOD-BY-LEETEN 2014/04/19:    def addblock_svd_update( U, s, V, A, force_orth ):
def addblock( U, s, V, A, force_orth ):
# MOD-BY-LEETEN 2014/04/19-END
    
    current_rank = U.shape[1];

    # P is an orthogonal basis of the column-space
    # of (I-UU')a, which is the component of "a" that is
    # orthogonal to U.
    m = np.dot(U.T, A);
    p = A - np.dot(U, m);

    # MOD-BY-LEETEN 2014/04/19:    P = la.orth( p );
    P = orth( p );
    # MOD-BY-LEETEN 2014/04/19-END
    # # MOD-BY-LEETEN 2014/04/20-FROM:
    # # p may not have full rank.  If not, P will be too small.  Pad
    # # with zeros.
    # P = np.pad(P, ((0,0), (0, p.shape[1] - P.shape[1])), 'constant', constant_values=(0, 0));
    # Ra = np.dot(P.T, p);
    # 
    # #  
    # # Diagonalize K, maintaining rank
    # #
    # 
    # z = np.zeros( m.shape );
    # K = np.vstack([np.hstack([np.diag(s), m]), np.hstack([z.T, Ra])]);
    # # MOD-BY-LEETEN 2014/04/20-TO:
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
    
    # # MOD-BY-LEETEN 2014/04/20-END
    # # 2 options to solve the SVD of k:
    # # Use SVD, which can be still storage consuming.    
    # [tUp, tsp, tVpT] = sparsela.svds( sparse.bsr_matrix(K), current_rank );
    # tUp = tUp[:, ::-1];
    # tsp = tsp[::-1];
    # tVp = tVpT[::-1, :].T;

    # Use PCA instead.
    # MOD-BY-LEETEN 2014/04/20:    tUp, tsp, tVp = regular_svd_by_pca(K, k = current_rank);
    tUp, tsp, tVp = regular_svd(K, k = current_rank);
    # MOD-BY-LEETEN 2014/04/20-END

    #
    # Now update our matrices!
    #

    sp = tsp;

    # this may not preserve orthogonality over
    # many repetitions.  See below.
    Up = np.dot(np.hstack([U, P]), tUp);   

    # Exploit structure to compute this fast: Vp = [ V Q ] * tVp;
    # MOD-BY-LEETEN 2014/04/21:    if( 0 != len(V) ): 
    if( 0 != V.size ):
    # MOD-BY-LEETEN 2014/04/21-END
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
        # # MOD-BY-LEETEN 2014/04/21-FROM:
        # UQ, UR = la.qr( Up, mode='economic' );
        # VQ, VR = la.qr( Vp, mode='economic' );
        # [tUp, tsp, tVp] = la.svd( np.dot(np.dot(UR, np.diag(sp)), VR.T) );
        # Up = np.dot(UQ, tUp);
        # Vp = np.dot(VQ, tVp);
        # sp = tsp;
        # # TO:
        Up, sp, Vp = re_orth(Up, sp, Vp);
        # # MOD-BY-LEETEN 2014/04/21-END

    return [Up, sp, Vp];

# TODO: Port the rank one update too.
# http://web.mit.edu/~wingated/www/scripts/rank_one_svd_update.m

# ADD-BY-LEETEN 2014/04/19-BEGIN

def my_addblock( U, s, V, A, force_orth ):
    new_n_cols = A.shape[1];
    V_padded = np.pad(V, ((0, new_n_cols), (0, 0)), 'constant', constant_values=(0, 0));
    B = np.pad(np.eye(new_n_cols, new_n_cols), ((V.shape[0], 0), (0, 0)), 'constant', constant_values=(0, 0));
    return update(U, s, V_padded, A, B, force_orth);
# ADD-BY-LEETEN 2014/04/19-END

# ADD-BY-LEETEN 2014/04/21-BEGIN
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
    UQ, UR = la.qr( Up, mode='economic' );
    VQ, VR = la.qr( Vp, mode='economic' );
    [tUp, tsp, tVp] = regular_svd( np.dot(np.dot(UR, np.diag(sp)), VR.T) );
    return np.dot(UQ, tUp), tsp, np.dot(VQ, tVp);
    
def addvector( U, s, V, a, force_orth ):
    V_padded = np.pad(V, ((0, 1), (0, 0)), 'constant', constant_values=(0, 0));
    b = np.pad(np.eye(1, 1), ((V.shape[0], 0), (0, 0)), 'constant', constant_values=(0, 0));
    return rank_one_update(U, s, V_padded, a, b, force_orth);    
# ADD-BY-LEETEN 2014/04/21-END
