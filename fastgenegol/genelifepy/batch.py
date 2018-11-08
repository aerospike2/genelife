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
# execute with "./activity.py ./genelifeAct.py"  

import sys
import cProfile
import numpy as np
import gmpy2 as gmp     # use gmpy2 package to allow efficient bit string operations (also for LEN>63)
from gmpy2 import mpz   # pip install gmpy2, doc see https://gmpy2.readthedocs.io/en/latest/mpz.html 
from copy import copy

def genelife_sub(args):
  """ returns value of simulation with these parameters (defn below)"""
  #%matplotlib notebook
  #%matplotlib inline

  p0 = float(args[0]); mutprob = float(args[1]); alpha = float(args[2]); LEN = float(args[3])
  initial1density = float(args[4]); NGC = float(args[5]); initmut = float(args[6]); neutral = float(args[7])
  LEN = int(LEN)
  NGC = int(NGC)
  
  # p0 = 0.1               # max prob of flip: compare with p0 = 0.0 to see advantage of model       (1)
  # mutprob = 0.1          # probability of single point mutation per replication                    (2)
  # alpha = 1.0            # exponential decay constant of flip prob with hamming distance           (3)
  # LEN = 63               # length of genome: LEN > 8 for current color display                     (4)
  # initial1density = 0.8  # initial density of ones in randomly set initial GoL pattern             (5)
  # NGC = 4                # no of initial gene centres in sequence space                            (6)
  # initmut = 0.2          # mutation prob for creating initial genes                                (7)
  # neutral = 1            # whether neutral or flip prob determined by select gene seq. (via nr 1s) (8)
  
  N = 128                  # size of array
  NG = 2^LEN -1            # max genome sequence as nr
  niterations = 400        # no of updates of grid in animation

  NC = LEN+1               # number of colors
  colormethod = 1          # 1 color by gene leading bits, 0 color by 1 + hamming(nbgenes)
# setup of color map:        black for 0, colors for 1 to LEN+1 or 257 for colormethod 0 or 1

  gradients = 0            # 1 add gradients in 2 key parameters, e.g. p0 in x and mutprob in y; 0 do not
  p0min = 0.0              # min value of p0 for gradient
  mutprobmin = 0.0         # minimum mutprob for gradient
  
  idact = {}               # initial list of ids for activity statistics
  iddone = {}

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
      """ gmp mutation of s """
      global totmut
      mcount = 0
      while np.random.random() < prob and mcount < LEN/2:
          mcount = mcount + 1
          pos = np.random.randint(0,LEN)
          if gmp.bit_test(s,pos): s = gmp.bit_clear(s,pos)
          else:                   s = gmp.bit_set(s,pos)
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
  
  
  def update(grid,genegrid,newgenegrid,idact,iddone):
  
      p0i = p0
      mutproby = mutprob
       
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
                  if total == 2 or total == 3: # Conway's rule unless genetic neighborhood says otherwise
                      nbgenes = one_neighbors(genegrid, nbs, i, j)
                      d = hamming(nbgenes) 
                      gs = rselect(nbgenes)                                   # genetic difference measure
                      if gradients: 
                          p0i = p0min + ((p0-p0min) * i) / float(N-1)
                          mutproby = mutprobmin+((mutprob-mutprobmin) * j) / float(N-1)
                      if neutral:
                          p = flipprob(d,p0i,alpha)
                      else:
                          n1av=float(LEN)/2.
                          p1 = p0i * max(0.,(gmp.popcount(gs)-n1av)/n1av)
                          p = flipprob(d,p1,alpha)                    
                      if total == 3:                    
                          if np.random.random() > p:                              # turn on if no flip
                              newgenegrid[i][j] = mutate(gs,mutproby)
                              newgrid[i, j] = 1
                      elif total == 2: # genetic neighborhood overrules Conway's rule (to keep off)
                          if np.random.random() < p:                              # turn on if flip
                              newgenegrid[i][j] = mutate(gs,mutproby)
                              newgrid[i, j] = 1
      # update data
      grid = newgrid
      iddone = {}
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
              iddone[x] = 1
              ## end activity computation 
              ####################################################
  
  
  def runvalue1(grid):
      """ density of 1s"""
      return grid.count(1)
  
  def runvalue2(grid, genegrid):
      """ overall or average local genetic diversity"""
      return 0
  
  def runvalue3(idact):
      """ integrated popoulation activity"""
      sum = 0
      for x in idact:
        if x != mpz(0):
          if x not in iddone:   # only live genes in last generation
            sum = sum + idact[x]
      return sum
  
  def runvalue4(grid):
      """ number of glider local patterns NYI"""
      return 0
  
  # populate grid with random integers and genes 
  seed = gmp.random_state(mpz(0x789abcdefedcba65))   # initialize seed for gmpy2 random operations
  # grid = np.random.randint(0, 2, N*N).reshape(N, N) # start with random grid of 0 or 1 values with equal probs
  grid = np.random.choice([0, 1], size=(N,N), p=[1.0-initial1density, initial1density]) # start with random grid of 0 or 1 values
  cgrid = grid.copy()    
  if not neutral:
      startgenecentres = [gmp.mpz_urandomb(seed, LEN) for i in range(NGC)]
      genegrid = [[(mutate(startgenecentres[np.random.choice([0, NGC-1])],initmut) if grid[i,j] else mpz(0)) for i in range(N)] for j in range(N)]
  else:
      genegrid = [[(gmp.mpz_urandomb(seed, LEN) if grid[i,j] else mpz(0)) for i in range(N)] for j in range(N)]
  newgenegrid = [[genegrid[i][j] for i in range(N)] for j in range(N)]
  
  for i in range(niterations):
    update(grid,genegrid,newgenegrid,idact,iddone)
  
  args = [float(aa) for aa in args]
  return args + [runvalue3(idact)]
  
if __name__ == '__main__':
  if len(sys.argv)==2:
    args = sys.argv[1]
    args = args.replace(","," ")
    args = args.replace("\\"," ")
    args = args.split()
  elif len(sys.argv)==9:
    args = sys.argv[1:]
  else:
    sys.exit("Usage:  batch.py arg1 ... arg8")
  foo = genelife_sub(args)
  for x in foo:
    print x,
  print
