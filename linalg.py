import numpy as np
import galois


GF2 = galois.GF(2)


def find_linear_dependence(vectors):
    '''Return a basis for the nullspace of the matrix formed by combining
    the input column vectors (over GF(2)).
    '''
    # Put the row vectors into a matrix
    # We want to find a solution to A^T v = 0
    # We do this by column operations on A^T or row operations on A
    A = np.array(vectors)
    ncols = len(vectors[0])
    nrows = len(vectors)
    # Create the augmented matrix A | I
    # I is an identity matrix with the same number of rows as A.
    # ones on the diagonal and zeroes elsewhere
    AI = np.concatenate((A, np.eye(nrows, dtype='int')), axis=1)
    AI_GF2 = GF2(AI)
    # Row reduce the A half of A | I
    AI_rr = AI_GF2.row_reduce(ncols=ncols)
    # Find the zero rows and put the corresponding basis vectors into an array
    basis = [AI_rr[i, ncols:].tolist() for i in range(nrows)
             if np.all(AI_rr[i, :ncols] == 0)]
    return basis
