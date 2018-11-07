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
libcd.genelife_update.argtypes = [ c_int, c_int, c_int]
libcd.initialize.restype = None
libcd.initialize.argtypes = [int_array, c_int, int_array, c_int]
libcd.initialize_planes.restype = None
libcd.initialize_planes.argtypes = [int_array, c_int]
libcd.countspecies1.restype = None
libcd.countspecies1.argtypes = [uint64_array, uint64_array, c_int]
libcd.countspecies.restype = None
libcd.countspecies.argtypes = None
libcd.countspecieshash.restype = None
libcd.countspecieshash.argtypes = None
#libcd.print_gol.restype = None
#libcd.print_gol.argtypes = [uint64_array, c_int, c_int]
#libcd.printscreen.restype = None
#libcd.printscreen.argtypes = [uint64_array, uint64_array, c_int, c_int]
libcd.get_histo.restype = None
libcd.get_histo.argtypes = [int_array, c_int]
libcd.get_curgol.restype = None
libcd.get_curgol.argtypes = [uint64_array,c_int]
libcd.get_curgolg.restype = None
libcd.get_curgolgstats.argtypes = [uint64_array,c_int]
libcd.get_curgolgstats.restype = None
libcd.get_curgolg.argtypes = [uint64_array,c_int]
libcd.get_stats.restype = None
libcd.get_stats.argtypes = [int_array, int_array, int_array, int_array, c_int]
libcd.get_activities.restype = None
libcd.get_activities.argtypes = [uint64_array, int_array, int_array]
libcd.colorgenes1.restype = None
libcd.colorgenes1.argtypes = [uint64_array, uint64_array, uint64_array, int_array, c_int]
libcd.colorgenes.restype = None
libcd.colorgenes.argtypes = [int_array, c_int]
libcd.set_colorfunction.restype = None
libcd.set_colorfunction.argtypes = [c_int]
libcd.setget_act_ymax.restype = c_int
libcd.setget_act_ymax.argtypes = [c_int]


def genelife_update(nsteps, nhist, nstat):
    return libcd.genelife_update(nsteps, nhist, nstat)

def get_curgol(gol):
    return libcd.get_curgol(gol, int(len(gol)))

def get_curgolg(golg):
    return libcd.get_curgolg(golg, int(len(golg)))

def get_curgolgstats(golgstats):
    return libcd.get_curgolgstats(golgstats, int(len(golgstats)))

def initialize(runparams,simparams):
    return libcd.initialize(runparams, len(runparams), simparams, len(simparams))

def initialize_planes(offsets):
    return libcd.initialize_planes(offsets, int(len(offsets)))

def countspecies1(gol, golg):
    return libcd.countspecies1(gol, golg,  len(golg))

def countspecies():
    return libcd.countspecies()

def countspecieshash():
    return libcd.countspecieshash()

#def print_gol( gol, N):
#    return libcd.print_gol( gol, N, len(gol))

#def printscreen( gol, golg, N):
#    return libcd.printscreen( gol, golg, N, len(gol))

def get_histo(gol):
    return libcd.getconfigs( histo, len(histo))

def get_stats(livesites,genotypes,stepstats,configstats,nstats):
    return libcd.get_stats(livesites,genotypes,stepstats,configstats,nstats)

def get_activities(actgenes,activities,ngenesp):
    return libcd.get_activities(actgenes,activities,ngenesp)

def colorgenes1(gol, golg, golgstats, cgolg):
    return libcd.colorgenes1( gol, golg, golgstats, cgolg, len(gol))

def colorgenes(cgolg):
    return libcd.colorgenes( cgolg, len(cgolg))

def set_colorfunction(colorfunctionval):
    return libcd.set_colorfunction(colorfunctionval)

def setget_act_ymax(ymax):
    return libcd.setget_act_ymax(ymax)
