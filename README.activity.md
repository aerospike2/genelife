Evolutionary-activity
=====================

Evolutionary activity is a visual and statistical approach to analyzing evolutionary processes.

Evolutionary activity is used to analyze an evolutionary process consisting of a population of entities that change over time.  Typically, the entities interact with each other, and as a result of the interactions, some entities may survive preferentially over others.  If entities have survival that is significantly longer than a randomly changing population with no selection, the entities acquire evolutionary activity.

$a_i^t = \Sum_{t'=0}^{t} c_i^t'$


Usage:

Configure your evolutionary simulation to print a line to standard out every time step.
This line should contain a unique identifier (string) for each type in the population,
and a number representing that type's strength in the population (simplest case: population
count of that type).  

Then run activity.py with arguments = your program with its arguments.
For example, if your program is evo, with arguments arg1 arg2 arg3, you would
visualize its evolutionary activity by running:

% activity.py evo arg1 arg2 arg3

For example, three timesteps of evo might look like:

A 1 a 4 aa 4 B 13 AA 3 ba 2 bb 1 aab 1 BB 79 BBAA 1 b 6 ab 2 AB 14 BA 42 
A 1 a 1 aa 8 B 5 AA 2 ba 1 bb 7 aab 1 BB 136 BBAA 1 b 3 BAA 3 ab 12 AB 3 BA 21 
A 9 a 1 aa 65 B 1 AA 13 ba 1 bb 17 aab 2 BB 9 BBAA 1 BBB 1 b 13 BAA 14 ab 8 AB 5 BA 1 

One gotcha: your program (evo in the example above) must flush stdout after writing the line each time step.

Example:

Run Nicholas Guttenberg's predator prey model and see its activity:

% activity.py guttt.py


To do:

2.  eliminate flat line clutter

3.  add graph windows for  total activity vs. time and new activity vs. time.

