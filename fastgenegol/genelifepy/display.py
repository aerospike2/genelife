#!/usr/local/bin/python

# disp_simulator is assumed to give output on stdout, one line per timestep,
# gene1 x1 y1 gene2 x2 y2 etc.
# for all live cells in the population that generation.
# if disp_simulator has args sim_arg1, sim_arg2, ...
# run with the following:
# display.py disp_simulator sim_arg1 sim_arg2 ...

# currently hardwired to take data from a 256x256 lattice, display to a 512x512 window (zoom x2)


from pygame import *
from math import *
from random import *
import sys
from time import sleep
import collections
import subprocess as sp
import os
from numpy import fromfile
from numpy import sin,pi
from os import system

Width = 512
Height = 512
#Width = 1024
#Height = 1024

# set window position
os.environ['SDL_VIDEO_WINDOW_POS'] = "%d,%d" % (20,20)

# rainbow colormap:
norm = lambda x: min(max(int((x+1)*128),0),255)
s = lambda t: sin(2*pi*t)
spec = lambda t: (norm(s(t*0.9+0.2)), norm(s(t*0.9+0.9)), norm(s(t*0.9+0.5)))
palette = tuple(spec(x/256.) for x in range(256))

def main():
    global ymax

    argv = sys.argv
    argv = argv[1:]
    sys.stderr.writelines(str(argv)+'\n')
    proc = sp.Popen(argv,
                    stdout = sp.PIPE)
    
    screen = display.set_mode([Width, Height])
    scr = surface.Surface([256,256], 0)
    display.set_caption("Gene Life")
    draw.rect(screen, [10, 10, 10],(0, 0 , Width, Height + 1), 0)
    display.update()

    gennum = 0
    mainloop = True
    while mainloop:
        for ev in event.get():
            # User presses QUIT-button.
            if ev.type == QUIT:
                mainloop = False 
            elif ev.type == KEYDOWN:
                # User presses ESCAPE-Key
                if ev.key == K_ESCAPE:
                    mainloop = False
        draw.rect(scr, [10, 10, 10],(0, 0 , Width/2, Height/2 + 1), 0)
        dat = proc.stdout.readline()
        dat = dat.split()
        # dat has (x,y,string) triples on each line.
        ldat = len(dat)
        if ldat == 0:
            return
        if ldat % 3:
            print "length of disp output not multiple of 3!"
            print '++++++++++ '
            print dat
            print '---------------\n'
        fofo = range(0,ldat,3)          # index into data to graph: fofo[i] of beginning of ith data chunk

        colors = {}
        mycol = 0
        for i in range(len(fofo)):
            mykey = dat[fofo[i]]
            X = int(dat[fofo[i]+1])
            Y = int(dat[fofo[i]+2])
            if mykey not in colors:
                colors[mykey] = palette[mycol%255]
                mycol += 29     # leap through rainbow
            scr.set_at((X,Y),colors[mykey])
        transform.scale2x(scr,screen)
        display.update()
        gennum = gennum +1

if __name__=='__main__':
    main()
    while True:
        for ee in event.get():
            if ee.type == QUIT:
                sys.exit()
            elif ee.type == KEYDOWN:
                if ee.key == K_ESCAPE or ee.key == K_q:
                    sys.exit()
                    
