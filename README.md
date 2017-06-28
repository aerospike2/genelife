# genelife
Genetic extension to Conway's Game of Life
by John McCaskill and Norman Packard

Description:

A simple Python/matplotlib implementation of Conway's Game of Life is
extended to include the influence of genes proliferating as directed by the game
and influencing the random innovations in the game, allowing long term interesting activity in GoL.

Individuals (with a gene) are associated only with \"live\" or \"on\" or \"1\" sites 
Individuals die and genes destroyed when a site \"dies\" i.e. is set to \"empty\" or \"off\" or \"0\"
Individuals replicate with point mutation in current version, later recombination
Normal replication is possible only to central empty (0) site and when 3 individuals are in neighborhood
The parent for mutation during replication is chosen randomly from neighbors
Conway's game of life rule is overridden stochastically only for empty sites with 2 or 3 neighbors present (on)
  The probability of rule override p=p0*e^-d decays exponentially with increasing hamming distance d of neighbors
  (i) In the case of three live neighbors, the central site remains dead and no replication happens
  (ii)In the case of two live neightbors, the central site comes alive and replication happens.

For p0 == 0, the occupied cells follow exactly Conway's game of life
  and the genes execute neutral selection from an initially random population
For p0 >0, the probability of departures from Conway's rules are greatest with monoclonal neighbors
  and become negligible if neighbors are distantly related
In this way, a feedback is created between pattern stagnation and innovation

With the current parameters, no degeneration to a set of non-communicating local structures occurs
The model is likely to be more interesting still with recombination than point mutation NYI

Starting point for python code was the modification of electronut.in by Takashi Ikegami.
