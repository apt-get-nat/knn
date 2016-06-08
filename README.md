# knn
A c implementation of the kdtree nearest neighbor search algorithm, for use as a DLL.

There are three methods meant to be called externally:

int SETUP( int dim, int N, double[] X)
    dim is the dimension of the problem (usually 2 or 3)
    N is the number of data locations
    X is a dim-wrapped array of data locations.
  This algorithm does not permute the input, but will create a KDTree in system memory.

ARRAYSEARCH(int dim, int N, double[] Y, int M, int[] K, double[] D)
    dim is the dimension of the problem -- the same as was given in SETUP
    N is the number of locations to find nearest neighbors of
    Y is the dim-wrapped array of N locations.
    M is the number of nearest neighbors to each requested
    K and D are pre-allocated arrays of length N*M.
  K and D will be permuted by the function. Upon return, K will contain the indices
  corresponding to nearest neighbors pulled from the X array in SETUP, and D will contain
  the corresponding cartesian distances.

INSPECT(int dim)
    dim is the dimension passed in SETUP
  Provides a (mostly) human-readable printout of the KDTree currently in memory.
