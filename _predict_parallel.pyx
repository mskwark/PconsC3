cimport numpy as np
cimport cython
from cython cimport parallel
cimport cython


@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline double predict_tree(long[:, ::1] trunk,  double[::1] leaf, double[::1] compare, double[::1] q) nogil:
    cdef double r
    cdef int i = 0, j

    while i >= 0:
        if q[trunk[0, i]] <= compare[i]:
            j = trunk[1, i]
        else:
            j = trunk[2, i]
        if j < 0:
            r = leaf[i]
        i = j
    return r


@cython.wraparound(False)
@cython.boundscheck(False)
def predict(long[:, :, ::1] trunks, double[:, ::1] leafs, double[:, ::1] compare, double[:, ::1] x,
            double[::1] predictions, int num_threads):
    if num_threads == 1:
        _predict_serial(trunks, leafs, compare, x, predictions)
    else:
        _predict_parallel(trunks, leafs, compare, x, predictions, num_threads)

    return predictions


@cython.wraparound(False)
@cython.boundscheck(False)
cdef void _predict_serial(long[:, :, ::1] trunks, double[:, ::1] leafs, double[:, ::1] compare,
                          double[:, ::1] x, double[::1] predictions) nogil:
    cdef int i, n, j, n_trees
    cdef double s
    n = x.shape[0]
    n_trees = trunks.shape[0]

    for i in xrange(n):
        s = 0.
        for j in xrange(n_trees):
            s += predict_tree(trunks[j, :, :], leafs[j, :], compare[j, :], x[i, :])
        predictions[i] = s / n_trees


@cython.wraparound(False)
@cython.boundscheck(False)
cdef void _predict_parallel(long[:, :, ::1] trunks, double[:, ::1] leafs, double[:, ::1] compare,
                            double[:, ::1] x, double[::1] predictions, int num_threads) nogil:
    cdef int i, n, j, n_trees
    cdef double s
    n = x.shape[0]
    n_trees = trunks.shape[0]

    for i in parallel.prange(n, num_threads=num_threads, chunksize=100, schedule='static', nogil=True):
        s = 0.
        for j in xrange(n_trees):
            s = s + predict_tree(trunks[j, :, :], leafs[j, :], compare[j, :], x[i, :])
        predictions[i] = s / n_trees
