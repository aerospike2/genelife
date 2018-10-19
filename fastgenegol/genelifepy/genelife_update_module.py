""" 
Wrapping a C library function that does update of long unsigned int arrays gol, golg
    input using the numpy.ctypeslib.
    Method gleaned from 
    http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html
"""

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
libcd.genelife_update.argtypes = [ c_int, c_int]
libcd.initialize.restype = None
libcd.initialize.argtypes = [int_array, c_int, int_array, c_int]
libcd.initialize_planes.restype = None
libcd.initialize_planes.argtypes = [int_array, c_int]
libcd.countspecies.restype = None
libcd.countspecies.argtypes = [uint64_array, uint64_array, int_array]
libcd.print_gol.restype = None
libcd.print_gol.argtypes = [uint64_array, c_int, c_int]
libcd.printscreen.restype = None
libcd.printscreen.argtypes = [uint64_array, uint64_array, c_int, c_int]
libcd.get_histo.restype = None
libcd.get_histo.argtypes = [uint64_array, c_int]
libcd.init_histo.restype = None
libcd.init_histo.argtypes = None
libcd.get_curgol.restype = None
libcd.get_curgol.argtypes = [uint64_array,c_int]
libcd.get_curgolg.restype = None
libcd.get_curgolg.argtypes = [uint64_array,c_int]
libcd.get_stats.restype = None
libcd.get_stats.argtypes = [int_array, int_array, c_int]
libcd.colorgenes.argtypes = [uint64_array, uint64_array, int_array, c_int]

def genelife_update(nsteps, histoflag):
    return libcd.genelife_update(nsteps, histoflag )

def get_curgol(gol):
    return libcd.get_curgol(gol, int(len(gol)))

def get_curgolg(golg):
    return libcd.get_curgolg(golg, int(len(golg)))

def initialize(runparams,simparams):
    return libcd.initialize(runparams, len(runparams), simparams, len(simparams))

def initialize_planes(offsets):
    return libcd.initialize_planes(offsets, len(offsets))

def countspecies(gol, golg, runparams):
    return libcd.countspecies(gol, golg,  runparams)

def print_gol( gol, N):
    return libcd.print_gol( gol, N, len(gol))

def printscreen( gol, golg, N):
    return libcd.printscreen( gol, golg, N, len(gol))

def get_histo(gol):
    return libcd.getconfigs( histo, len(histo))

def get_stats(livesites,genotypes,nstats):
    return libcd.get_stats(livesites,genotypes,nstats)

def colorgenes(gol, golg, cgolg):
    return libcd.colorgenes( gol, golg, cgolg, len(gol))
