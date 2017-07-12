
<script type="text/javascript" src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_HTML"></script>

# genelife

Genetic extension to Conway's Game of Life
by John McCaskill and Norman Packard

## Introduction

Conway's Game of Life is extended to include the influence of genes
proliferating as part of the game dynamics.  The genes cause a
probabalistic departure from the Game of Life only for certain
neighborhood configurations.  The resulting Genetic Game of Life has
random innovations as a result of the genetic dynamics, enabling long
term interesting dynamics, in contrast to the simple fixed point and
oscillating local states characteristic of the Game of Life.

## Conway's Game of Life

Conway's Game of Life is a two dimensional cellular automaton rule.
States are configurations of ones and zeros on a two dimensional
lattice.  States change in time according to a local rule applied
simultaneously to all sites.  The local rule for a given site is a
function of that site's value as well as the value of the eight nearest
neighboring sites.  Using the convention that sites with value one are
considered live, and sites with value zero are considered dead, the
local rule may be stated as follows:

1. Any live cell with fewer than two live neighbours dies, as if
caused by underpopulation.

2. Any live cell with two or three live neighbours lives on to the
next generation.

3. Any live cell with more than three live neighbours dies, as if by
overpopulation.

4. Any  dead cell with  exactly three  live neighbours becomes  a live
cell, as if by reproduction.

## The Genetic Game of Life.

Live sites are considered as an individual whose identity is
determined by an associated genome consisting of a bit string.
Individuals die and its genes are destroyed when the site dies.  When
individuals replicate in step 4 of the Conway's rule, a parent is
chosen randomly from the 3 live neighbors, and the genes copied from
that parent to the new individuals genome.  The genome may change,
with point mutation in current version (recombination to come in the
future).

In order to give a feedback between genes and the CA time development, 
Conway's game of life rule is overridden stochastically only for empty
sites with 2 or 3 neighbors present (on):
(i) In the case of three live neighbors, the central site remains dead
and no replication happens
(ii) In the case of two live neightbors, the central site comes alive
and replication happens.

The probability of rule override $p=p0*e^-alphha*d$ decays exponentially with
increasing hamming distance d of neighbors.
For p0 == 0, the occupied cells follow exactly Conway's game of life
and the genes execute neutral selection from an initially random
population

For p0 >0, the probability of departures from Conway's rules are
greatest with monoclonal neighbors and become negligible if neighbors
are distantly related.

In this way, a feedback is created between pattern stagnation and innovation.

With the current parameters, no degeneration to a set of
non-communicating local structures occurs. The model is likely to be
more interesting still with recombination than point mutation NYI

## Versions

Starting point for python code was the modification of electronut.in
by Takashi Ikegami.

* genelife.ipynb:  ipython notebook.
* genelife.py: version run at the command line with "python genelife.py"

### Evolutionary Activity

Activity statistics assign an activity counter to each genome.  At
each time step, the counter is incremented for every occurance of the
genenome in the population over the entire lattice.  In other words,
the counter is incremented by the genome population for that time
step.  As the genome proliferates through reproduction, its activity
increases.  When the genome is mutated during a reproductive event,
the new genome has its own counter, initialized to zero (unless the
genome already exists, in which case its counter is whatever it is).

* genelifeAct.py: version that outputs activity statistics, run with
  "activity.py genelifeAct.py" (cf. README.activity.md)

### Batch runs for optimization

A version for batch exploration and optimization exposes parameters with a 
command line argument.  

The arguments are (with default values listed):

1.  p0 = 0.1:  max prob of flip: compare with p0 = 0.0 to see advantage of model       
2.  mutprob = 0.1:  probability of single point mutation per replication                    
3.  alpha = 1.0:  exponential decay constant of flip prob with hamming distance           
4.  LEN = 63:  length of genome: LEN > 8 for current color display                     
5.  initial1density = 0.8:  initial density of ones in randomly set initial GoL pattern             
6.  NGC = 4:  no of initial gene centres                                              
7.  initmut = 0.2:  mutation prob for creating initial genes                                
8.  neutral = 1:  whether neutral or with p0 determined by select gene seq. (via nr 1s)   

An example call using these arguments is:

```
batch.py 0.5 0.5 0 15 0.7 8 0.125 0
```

For reasons of interfacing to GNU parallel, it is also convenient to
enable passing arguments as a single string, .e.g:

```
batch.py "0.5 0.5 0 15 0.7 8 0.125 0"
```

## Optimization workflow with PDT

* set up experimental space, e.g. with expdef.csv
```
Name,Value.1,Value.2,Value.3,Value.4,Value.5,Value.6,Value.7,Value.8,Value.9,Value.10,Value.11
p0,0,0.001,0.002,0.004,0.008,0.016,0.031,0.062,0.125,0.25,0.5
mutprob,0,0.001,0.002,0.004,0.008,0.016,0.031,0.062,0.125,0.25,0.5
alpha,0,0.5,1,1.5,2,2.5,3,,,,
LEN,7,15,31,63,127,,,,,,
initial1density,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,,
NGC,1,2,3,4,8,16,32,,,,
initmut,0,0.008,0.016,0.031,0.062,0.125,0.25,0.5,,,
neutral,0,1,,,,,,,,,
```
* run PDT at protolife.com
* when a design is produced, copy and paste into a file, e.g. `design-001`
* run all the designs with parallel:
```
cat design-001 | parallel batch.py {} > response-001.tmp
```
NB:  you can ssh to a compute server that has lots of cores to run
the designs faster.
* process the output of parallel to match response calculation with
  correct designs:
```
process.py design.001 response-001.tmp > response-001
```

* paste the responses back into the PDT web page table.NB: to get just
  the final response column you could say
```
cat response-001 | awk '{print $NF}'
```
* Hit the button to produce the next generation of designs.
* repeat.