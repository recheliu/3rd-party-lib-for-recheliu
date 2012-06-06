Downloaded from: http://www.mit.edu/~wingated/resources.html

Quoted from the original website:

Incremental SVD updates

    Given the thin SVD of a matrix (X=USV'), update it in a number of interesting ways, while preserving the rank of the result.
    svd_update.m - update the SVD to be [X + A'*B]=Up*Sp*Vp' (a general matrix update).
    addblock_svd_update.m - update the SVD to be [X A] = Up*Sp*Vp' (add columns [ie, new data points])
    rank_one_svd_update.m - update the SVD to be [X + a*b'] = Up*Sp*Vp' (that is, a general rank-one update. This can be used to add columns, zero columns, change columns, recenter the matrix, etc. ).
     
    See Matthew Brand, "Fast low-rank modifications of the thin singular value decomposition". 


Teng-Yok Lee
3:44 PM 5/18/2012