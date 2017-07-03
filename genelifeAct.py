#!/usr/bin/env python

################################################################################
# genelife.py
#
# by John S. McCaskill and Norman H. Packard
#
# Based on Conway game of Life code by electronut.in as modified by Takashi Ikegami 

# Description:
#
# A simple Python/matplotlib implementation of Conway's Game of Life is
# extended to include the influence of genes proliferating as directed by the game
# and influencing the random innovations in the game
#
# 1. Individuals (with a gene) are associated only with "live" or "on" or "1" sites 
# 2. Individuals die and genes destroyed when a site "dies" i.e. is set to "empty" or "off" or "0"
# 3. Individuals replicate with point mutation in current version, later with recombination
# Normal replication is possible only to a central empty (0) site and when 3 individuals are in neighborhood
# The parent for mutation during replication is chosen randomly from neighbors
# 4. Conway's game of life rule is overridden stochastically only for empty sites with 2 or 3 neighbors present (on)
#   The probability of rule override p=p0*e^(-a*d) decays exponentially with increasing hamming distance d of neighbors
#   (i) In the case of three live neighbors, the central site remains dead and no replication happens
#   (ii)In the case of two live neightbors, the central site comes alive and replication happens.
#
# For p0 == 0, the occupied cells follow exactly Conway's game of life
#   and the genes execute neutral selection from an initially random population
# For p0 >0, the probability of departures from Conway's rules are greatest with monoclonal neighbors
#   and become negligible if neighbors are distantly related
# In this way, a feedback is created between pattern stagnation and innovation, 
#   and between replication and mutational diversity.
# With the current parameters, no degeneration to a set of non-communicating local structures occurs
# The model is likely to be more interesting still with recombination than point mutation NYI
################################################################################

import sys
import cProfile
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpcolors
import matplotlib.animation as animation
from matplotlib.colors import ListedColormap
import matplotlib
import gmpy2 as gmp     # use gmpy2 package to allow efficient bit string operations (also for LEN>63)
from gmpy2 import mpz   # pip install gmpy2, doc see https://gmpy2.readthedocs.io/en/latest/mpz.html 
from copy import copy
#%matplotlib notebook
#%matplotlib inline
N = 128         # size of array
LEN = 63        # length of genome: LEN > 8 for current color display
NG = 2^LEN -1   # max genome sequence as nr
NC = LEN+1      # number of colors
p0 = 0.08        # max prob of flip: compare with p0 = 0.0 to see advantage of model
alpha = 1       # exponential decay constant of flip prob with hamming distance
mutprob = 0.3   # probability of single point mutation per replication
colormethod = 1 # 1 color by gene leading bits, 0 color by 1 + hamming(nbgenes)

# setup of color map : black for 0, colors for 1 to LEN+1 or 257 for colormethod 0 or 1
#-----------------------------------------------------------------------------------------------------------
def rand_cmap(nlabels, type='bright', first_color_black=True, last_color_black=False, verbose=True):
    """
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks
    :param nlabels: Number of labels (size of colormap)
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :param first_color_black: Option to use first color as black, True or False
    :param last_color_black: Option to use last color as black, True or False
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :return: colormap for matplotlib
    Author of this color function: Delestro, stackoverflow or https://github.com/delestro/rand_cmap
    """
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np

    if type not in ('bright', 'soft'):
        print ('Please choose "bright" or "soft" for type')
        return
    if verbose:
        print('Number of labels: ' + str(nlabels))

    # Generate color map for bright colors, based on hsv
    if type == 'bright':
        randHSVcolors = [(np.random.uniform(low=0.0, high=1),
                      np.random.uniform(low=0.2, high=1),
                      np.random.uniform(low=0.9, high=1)) for i in xrange(nlabels)]    
        randRGBcolors = []  
        for HSVcolor in randHSVcolors:  # Convert HSV list to RGB
            randRGBcolors.append(colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2]))
        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]
        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Generate soft pastel colors, by limiting the RGB spectrum
    if type == 'soft':
        low = 0.6
        high = 0.95
        randRGBcolors = [(np.random.uniform(low=low, high=high),
                      np.random.uniform(low=low, high=high),
                      np.random.uniform(low=low, high=high)) for i in xrange(nlabels)]
        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]
        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Display colorbar
    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))
        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)
        cb = colorbar.ColorbarBase(ax, cmap=random_colormap, norm=norm, spacing='proportional', ticks=None,
                               boundaries=bounds, format='%1i', orientation=u'horizontal')

    return random_colormap
#-----------------------------------------------------------------------------------------------------------
if colormethod == 1:
    my_cmap = rand_cmap(257, type='bright', first_color_black=True, last_color_black=False, verbose=False)
else:
    my_cmap = matplotlib.cm.get_cmap('rainbow')  #was brg or spring or rainbow or PuBuGn
    my_cmap.set_under('black')  # use with vmin = 0.001 to set only the 0 integer state to black


def hamming(slist):
    """ extends hamming distance to many sequences using gmpy2"""
    l = len(slist)
    ds = 0
    for i in range(l-2):
           ds = ds + gmp.hamdist(slist[i],slist[i+1])
    if l>1: return int(ds//(l-1))
    else: return 0
           
def flipprob(d,p0,alpha):
    """flip probability calculated with exponential decay (const alpha)
       with hamming distance from max. value p0"""
    return p0*np.exp(-alpha*d)

def rselect(slist):
    """ select random element of list slist"""
    return slist[np.random.randint(0,len(slist)-1)]
    
def mutate(s,prob):
    """ gmp mutation of s: currently at most one point mutation,
        because only considering low rates """
    if np.random.random() < prob:
        pos = np.random.randint(0,LEN)
        if gmp.bit_test(s,pos): gmp.bit_clear(s,pos)
        else:                   gmp.bit_set(s,pos)
    return s

def neighbors_np(g, i, j):
    """ make list of 8 g neighbors of a site (i,j) with periodic bcs
    for numpy 2d array """
    if j > 1 and j < N-1:
        if i > 1 and i < N-1:   # separate calculation without modulo for speed
            return  [ g[i  , j-1], g[i  , j+1], 
                      g[i-1, j  ], g[i+1, j  ], 
                      g[i-1, j-1], g[i-1, j+1], 
                      g[i+1, j-1], g[i+1, j+1]]
        else:                   # modulo in i needed
            return  [ g[i,    j-1], g[i,       j+1], 
                      g[(i-1)%N, j  ], g[(i+1)%N, j  ], 
                      g[(i-1)%N, j-1], g[(i-1)%N, j+1], 
                      g[(i+1)%N, j-1], g[(i+1)%N, j+1]]
    else:                       # general case
        return  [ g[i,      (j-1)%N], g[i      ,(j+1)%N], 
                  g[(i-1)%N, j     ], g[(i+1)%N,j      ], 
                  g[(i-1)%N,(j-1)%N], g[(i-1)%N,(j+1)%N], 
                  g[(i+1)%N,(j-1)%N], g[(i+1)%N,(j+1)%N]]

def one_neighbors(gg, nbs, i, j):
    """" make list of gg neighbors of a site (i,j) at which grid value stored in nbs is True
        and for 1d array of 8 neighbors nbs"""
    onenbs = []
    k=0
    if j > 1 and j < N-1 and i > 1 and i < N-1: # separate calculation without modulo for speed
        if nbs[k]: onenbs.append(gg[i][j-1])
        k=k+1
        if nbs[k]: onenbs.append(gg[i][j+1])
        k=k+1
        if nbs[k]: onenbs.append(gg[i-1][j])
        k=k+1
        if nbs[k]: onenbs.append(gg[i+1][j] )
        k=k+1
        if nbs[k]: onenbs.append(gg[i-1][j-1])
        k=k+1
        if nbs[k]: onenbs.append(gg[i-1][j+1])
        k=k+1
        if nbs[k]: onenbs.append(gg[i+1][j-1])
        k=k+1
        if nbs[k]: onenbs.append(gg[i+1][j+1])
    else:                                       # general case with modulo
        if nbs[k]: onenbs.append(gg[i][(j-1)%N])
        k=k+1
        if nbs[k]: onenbs.append(gg[i][(j+1)%N])
        k=k+1
        if nbs[k]: onenbs.append(gg[(i-1)%N][j])
        k=k+1
        if nbs[k]: onenbs.append(gg[(i+1)%N][j] )
        k=k+1
        if nbs[k]: onenbs.append(gg[(i-1)%N][(j-1)%N])
        k=k+1
        if nbs[k]: onenbs.append(gg[(i-1)%N][(j+1)%N])
        k=k+1
        if nbs[k]: onenbs.append(gg[(i+1)%N][(j-1)%N])
        k=k+1
        if nbs[k]: onenbs.append(gg[(i+1)%N][(j+1)%N])
    return onenbs

def colorgrid(N,LEN,colormethod):
    """ colors array according to grid and genegrid using colormethod"""
    global grid,cgrid,genegrid
    if colormethod: # color by gene types            
        for i in range(N):                 # make this simple array copy operation more efficient
            for j in range(N): 
                if grid[i,j]:
                    cgrid[i,j] = 1 + gmp.c_div_2exp(genegrid[i][j],LEN-8) 
                else:
                    cgrid[i,j] = 0 
    else:           # color by hamming distance of neighbor sites with grid one value
        for i in range(N):                 # make this simple array copy operation more efficient
            for j in range(N): 
                if grid[i,j]:     
                    nbs = neighbors_np(grid, i, j)
                    nbgenes = one_neighbors(genegrid,nbs,i,j)
                    cgrid[i,j] = 1 + hamming(nbgenes) 
                else:
                    cgrid[i,j] = 0 
    return

idact = {}

def update(data):
    global grid, cgrid
    global genegrid, newgenegrid
     
    # copy grid since we need to update grid in parallel using old neighbors
    newgrid = grid.copy()
    for i in range(N):                 # make this simple array copy operation more efficient
        for j in range(N):
            newgenegrid[i][j] = genegrid[i][j]

    for i in range(N):
        for j in range(N):
            # compute 8-neighbor sum using toroidal boundary conditions
            nbs = neighbors_np(grid, i, j)
            total = sum(nbs)
            # apply Conway's rules with stochastic override for off sites with 2 or 3 neighbors
            if grid[i, j]  == 1: # cell value is on, "live" site
                if (total < 2) or (total > 3): # value set to off (zero), gene "dies"
                    newgrid[i, j] = 0
                    newgenegrid[i][j] = mpz(0)
                # else:  total is two or three: value stays on (live), nothing to do, no replication
            else:              # cell value is off, "dead" site
                if total == 3: # Conway's rule unless genetic neighborhood says otherwise
                    nbgenes = one_neighbors(genegrid, nbs, i, j)
                    d = hamming(nbgenes)                                    # genetic difference measure
                    p = flipprob(d,p0,alpha)                                # probability of flip override
                    if np.random.random() > p:                              # turn on if no flip
                        newgenegrid[i][j] = mutate(rselect(nbgenes),mutprob)
                        newgrid[i, j] = 1

                elif total == 2: # genetic neighborhood overrules Conway's rule (to keep off)
                    nbgenes = one_neighbors(genegrid, nbs, i, j)
                    d = hamming(nbgenes)                                    # genetic difference measure
                    p = flipprob(d,p0,alpha)                                # probability of flip override
                    if np.random.random() < p:                              # turn on if flip
                        newgenegrid[i][j] = mutate(rselect(nbgenes),mutprob)
                        newgrid[i, j] = 1
    # update data
    grid = newgrid
    for i in range(N):                 # make this simple array copy operation more efficient
        for j in range(N):
            genegrid[i][j] = newgenegrid[i][j]
            ####################################################
            ## here is the activity computation 
            x = genegrid[i][j]
            if x in idact:
                idact[x] += 1
            else:
                idact[x] = 1

    strout = ''
    for x in idact:
        strout = strout + str(x) + ' ' + str(idact[x]) + ' '
    print strout
    sys.stdout.flush()
    msg = 'length(idact) = ' + str(len(idact)) + '\n'
    sys.stderr.write(msg)
    ## end activity computation 
    ####################################################

    colorgrid(N,LEN,colormethod)        
    mat.set_data(cgrid)
    # hist = np.histogram(cgrid, bins=257)
    # print hist
    return [mat]

# populate grid with random integers and genes 
seed = gmp.random_state(mpz(0x789abcdefedcba65))   # initialize seed for gmpy2 random operations
grid = np.random.randint(0, 2, N*N).reshape(N, N) # start with random grid of 0 or 1 values
cgrid = grid.copy()    
genegrid = [[(gmp.mpz_urandomb(seed, LEN) if grid[i,j] else mpz(0)) for i in range(N)] for j in range(N)]
colorgrid(N,LEN,colormethod)
newgenegrid = [[genegrid[i][j] for i in range(N)] for j in range(N)]

# set up animation
fig, ax = plt.subplots()
mat = ax.matshow(cgrid, cmap=my_cmap, vmin=0.01, vmax=257)  # was vmax = LEN+1
ani = animation.FuncAnimation(fig, update, interval=10,
                              save_count=1, frames=1000, repeat = False)

plt.show()