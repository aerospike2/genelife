""" Wrapping a C library function that does update of long unsigned int arrays gol, golg
    input using the numpy.ctypeslib. """

import numpy as np
import numpy.ctypeslib as npct
from ctypes import c_int

# input type for the genelife_update function
# must be a long unsigned int array, with single dimension that is contiguous
uint64_array = npct.ndpointer(dtype=np.uint64, ndim=1, flags='CONTIGUOUS')
int_array = npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')

# load the library, using numpy mechanisms
libcd = npct.load_library("libgenelife", ".")

# setup the return types and argument types
libcd.genelife_update.restype = None
libcd.genelife_update.argtypes = [uint64_array, uint64_array, c_int, c_int, int_array, c_int, c_int]
libcd.initialize.restype = None
libcd.initialize.argtypes = [uint64_array, int_array, c_int, c_int]
libcd.initialize_genes.restype = None
libcd.initialize_genes.argtypes = [uint64_array, uint64_array, int_array, c_int, c_int]
libcd.initialize_planes.restype = None
libcd.initialize_planes.argtypes = [uint64_array, int_array, c_int, c_int]
libcd.countspecies.restype = None
libcd.countspecies.argtypes = [uint64_array, int_array, c_int, c_int]
libcd.print_gol.restype = None
libcd.print_gol.argtypes = [uint64_array, c_int, c_int]
libcd.printscreen.restype = None
libcd.printscreen.argtypes = [uint64_array, uint64_array, c_int, c_int]

def genelife_update(gol, golg, log2N, ndsteps, simparams):
    return libcd.genelife_update(gol, golg, log2N, ndsteps, simparams, len(gol), len(simparams) )

def initialize(gol,simparams):
    return libcd.initialize(gol, simparams, len(gol),len(simparams))

def initialize_genes(golg, gol, simparams):
    return libcd.initialize_genes(golg, gol, simparams, len(golg), len(simparams))

def initialize_planes(planeparams):
    return libcd.initialize_planes(offsets, planeparams,len(offsets), len(planeparams))

def countspecies(golg, simparams):
    return libcd.countspecies(golg,  simparams, len(golg), len(simparams))

def print_gol( gol, N):
    return libcd.print_gol( gol, N, len(gol))

def printscreen( gol, golg, N):
    return libcd.printscreen( gol, golg, N, len(gol))
