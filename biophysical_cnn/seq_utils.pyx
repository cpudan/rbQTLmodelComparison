#cython: language_level=3

# %%cython
import numpy as np
np.get_include()
cimport cython 
cimport numpy as np 

REVERSER=bytes.maketrans(b'AGCT',b'TCGA')

def reverse_complement(seq):
    return seq.translate(REVERSER)[::-1]

cdef dict bases={ 'A':<int>0, 'C':<int>1, 'G':<int>2, 'T':<int>3 } 

@cython.boundscheck(False) 
def one_hot( str string ):
    cdef np.ndarray[np.float32_t, ndim=2] res = np.zeros( (4,len(string)), dtype=np.float32 )
    cdef int j
    for j in range(len(string)):
        if string[j] in bases: # bases can be 'N' signifying missing: this corresponds to all 0 in the encoding
            res[ bases[ string[j] ], j ]=float(1.0)
    return(res)