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

# A wrapper of svd.
# Otherwise, scipy.sparse.linals's svd cannot handle the case that the rank is equal to the min(#rows, #col).  
def regular_svd(A, k=0):
    if( k <= 0 or k >= min(A.shape) ):
        return la.svd(A, full_matrices=False);
    else:    
        return sparsela.svds(sparse.bsr_matrix(A), k = k);    

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

    P = la.orth( p );
    # Pad with zeros if p does not have full rank.
    P = np.hstack([P, np.zeros(P.shape[0], p.shape[1] - P.shape[1])]);

    Ra = np.dot(P.T, p);

    #: Q: An orthogonal basis of the column-space of (I-VV')b.
    n = np.dot(V.T, B);
    q = B - np.dot(V, n);

    Q = la.orth( q );
    # Pad with zeros if q does not have full rank.
    Q = np.hstack([Q, np.zeros(Q.shape[0], q.shape[1] - Q.shape[1])]); 

    Rb = np.dot(Q.T, q);
 
    # Diagonalize K, maintaining rank.
    z = np.zeros( m.shape );
    z2 = np.zeros( [m.shape[1], m.shape[1]] );

    K = np.vstack([np.hstack([np.diag(s), z]), np.hstack([z.T, z2])]) + np.vstack([m, Ra]) * np.hstack(n, Rb);

    # tUp, tsp, tVpT = sparsela.svds(sparse.bsr_matrix(K), k=current_rank );
    # tUp = tUp[:, ::-1];
    # tSp = np.diag(tsp[::-1]); 
    # tVp = tVpT[::-1, :].T;
    tUp, tsp, tVp = regular_svd_by_pca(K, k = current_rank);
  
    # update our matrices.

    sp = tsp;
    Up = np.dot(np.hstack([U, P]), tUp);
    Vp = np.dot(np.hstack([V, Q]), tVp);

    # The above rotations may not preserve orthogonality, so we explicitly
    # deal with that via a QR plus another SVD.  In a long loop, you may
    # want to force orthogonality every so often.

    if ( force_orth ):
        UQ, UR = la.qr( Up, model='economic' );
        VQ, VR = la.qr( Vp, model='economic' );
        [tUp, tsp, tVp] = la.svd( np.dot(np.dot(UR, np.diag(sp)), VR.T) );
        Up = np.dot(UQ, tUp);
        Vp = np.dot(VQ, tVp);
        sp = tsp;
  
    return [Up, sp, Vp];

#! Given X = USV', compute SVD when appending new columns A to X
#
# Based on http://web.mit.edu/~wingated/www/scripts/addblock_svd_update.m
#

def addblock_svd_update( U, s, V, A, force_orth ):
    current_rank = U.shape[1];

    # P is an orthogonal basis of the column-space
    # of (I-UU')a, which is the component of "a" that is
    # orthogonal to U.
    m = np.dot(U.T, A);
    p = A - np.dot(U, m);

    P = la.orth( p );
    # p may not have full rank.  If not, P will be too small.  Pad
    # with zeros.
    P = np.pad(P, ((0,0), (0, p.shape[1] - P.shape[1])), 'constant', constant_values=(0, 0));

    Ra = np.dot(P.T, p);

    #  
    # Diagonalize K, maintaining rank
    #

    z = np.zeros( m.shape );
    K = np.vstack([np.hstack([np.diag(s), m]), np.hstack([z.T, Ra])]);

    # % TEST-MOD-FROM: 
    
    # # 2 options to solve the SVD of k:
    # # Use SVD, which can be still storage consuming.    
    # [tUp, tsp, tVpT] = sparsela.svds( sparse.bsr_matrix(K), current_rank );
    # tUp = tUp[:, ::-1];
    # tsp = tsp[::-1];
    # tVp = tVpT[::-1, :].T;

    # Use PCA instead.
    tUp, tsp, tVp = regular_svd_by_pca(K, k = current_rank);

    #
    # Now update our matrices!
    #

    sp = tsp;

    # this may not preserve orthogonality over
    # many repetitions.  See below.
    Up = np.dot(np.hstack([U, P]), tUp);   

    # Exploit structure to compute this fast: Vp = [ V Q ] * tVp;
    if( 0 != len(V) ): 
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
        UQ, UR = la.qr( Up, mode='economic' );
        VQ, VR = la.qr( Vp, mode='economic' );
        [tUp, tsp, tVp] = la.svd( np.dot(np.dot(UR, np.diag(sp)), VR.T) );
        Up = np.dot(UQ, tUp);
        Vp = np.dot(VQ, tVp);
        sp = tsp;

    return [Up, sp, Vp];

# TODO: Port the rank one update too.
# http://web.mit.edu/~wingated/www/scripts/rank_one_svd_update.m
