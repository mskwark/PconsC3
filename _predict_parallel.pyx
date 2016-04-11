cimport numpy as np
cimport cython
from cython cimport parallel

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline double predict_tree(long[:, ::1] tree,  double[:, :, ::1] leaf, double[::1] q) nogil:
    cdef double r = 0
    cdef int i = 0, j, saved_i

    while i >= 0:
        if q[tree[0, i]] <= tree[1, i]:
            j = tree[2, i]
        else:
            j = tree[3, i]
        if j < 0:
            saved_i = i
        i = j
    r = leaf[saved_i, 0, 1] / (leaf[saved_i, 0, 1] + leaf[saved_i, 0, 0] + 0.0001)
    return r



def predict(long[:, :, ::1] trees, double[:, :, :, ::1] leafs, double[:, ::1] X,
            double[::1] predictions, int num_threads):
    if num_threads == 1:
        _predict_serial(trees, leafs, X, predictions)
    else:
        _predict_parallel(trees, leafs, X, predictions, num_threads)

    return predictions


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void _predict_serial(long[:, :, ::1] trees, double[:, :, :, ::1] leafs, double[:, ::1] X,
            double[::1] predictions) nogil:
    cdef int i, n, j, n_trees
    cdef double s
    n = X.shape[0]
    n_trees = trees.shape[0]

    for i in xrange(n):
        s = 0.
        for j in xrange(n_trees):
            s += predict_tree(trees[j, :, :], leafs[j, :, :, :], X[i, :])
        predictions[i] = s / n_trees


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void _predict_parallel(long[:, :, ::1] trees, double[:, :, :, ::1] leafs, double[:, ::1] X,
            double[::1] predictions, int num_threads) nogil:
    cdef int i, n, j, n_trees
    cdef double s
    n = X.shape[0]
    n_trees = trees.shape[0]

    for i in parallel.prange(n, num_threads=num_threads, chunksize=20, schedule='static'):
        s = 0.
        for j in xrange(n_trees):
            # s += means reduction variable, don't use it
            s = s + predict_tree(trees[j, :, :], leafs[j, :, :, :], X[i, :])
        predictions[i] = s / n_trees
