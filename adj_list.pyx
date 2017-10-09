cimport numpy as np
import cython

@cython.wraparound(False)
@cython.nonecheck(False)
@cython.boundscheck(False)
def isNeighbor2(np.ndarray s1, np.ndarray s2, size_t length):
    """
    Args:
        s1: 1D array of ints
        s2: 1D array of ints
        length: length of s1 and s2 (must be the same)
    
    Returns:
        True if s1 and s2 have exactly 1 difference
        False otherwise (0, >1 differences)
    """
    
    cdef int mismatches = 0
    cdef size_t i
    for i in range(length):
        if s1[i] == s2[i]:
            continue
        else:
            mismatches += 1
        if mismatches > 1:
            return False
    if mismatches == 1:
        return True
    else:
        return False
