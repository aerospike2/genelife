#!/usr/local/bin/python

# like graphactivity.py, but try and use subprocess module so that execution is simply
# activity.py simulator sim_arg1 sim_arg2 ...


from pygame import *
from math import *
from random import *
import sys
from time import sleep
import collections
import subprocess as sp

# import pdb;

# set to None if you want to skip writing data out...
datfile = '/tmp/activity.dat'
if datfile is not None:
    datout = open(datfile,'w')
    
maxact = 5000                            # max number of activity traces in graph...
Width = 1000
Height = 500
ymax = 1000
ncount = 0                          # number of calls to trace

def trace(screen, colvalvec):            # eachcolval = [key,activityvalue]
                                         #        self.traceinit()
    cnt = 0
    cnt1 = 0
    yvals = []
    cols = []
    global ymax
    global ncount
    #        pdb.set_trace()
    for colval in colvalvec:            # separate color from yvals
        cols.append(colval[0])
        yvals.append(colval[1])
        cnt += 1
        if cnt>maxact:
            print "too many activity points..."
            break
    # do the scroll:
    if ncount<Width:      # first, don't scroll
        yvals = [Height - y * Height / ymax
                 for y in yvals]
        cnt1 = 0
        for i in range(len(yvals)):
            x = ncount
            y = yvals[i]
            col = cols[i]
            screen.set_at((x,y),col)
    else:                                 #  then scroll when full
        yvals = [Height - y * Height / ymax
                 for y in yvals]
        screen.scroll(-1,0)    # -1 => 1 pixel to left
        draw.line(screen,[10,10,10], # grey line first
                  (Width-1,Height),(Width-1,0))
        for i in range(len(yvals)):
            x = Width-1
            y = yvals[i]
            col = cols[i]
            screen.set_at((x,y),col)
    ncount += 1
    display.update()


import os
from numpy import fromfile
from numpy import sin,pi
from os import system

# set window position
os.environ['SDL_VIDEO_WINDOW_POS'] = "%d,%d" % (20,20)

# rainbow colormap:
norm = lambda x: min(max(int((x+1)*128),0),255)
s = lambda t: sin(2*pi*t)
spec = lambda t: (norm(s(t*0.9+0.2)), norm(s(t*0.9+0.9)), norm(s(t*0.9+0.5)))
palette = tuple(spec(x/256.) for x in range(256))


def omain():
    system("rm -f /tmp/activity; mkfifo -m 666 /tmp/activity")
    pipename = "/tmp/activity"
    print 'opening',pipename
    pipefd = open(pipename,"r")
    print 'opened',pipename
    while True:
        dat = pipefd.readline()
        dat = [int(x) for x in dat.split()]
        print dat
        print '\n---------------\n'


import pygame
def ooomain():

#     # foo = Graph(0,250,Width,Height,"Activity")
#     screen = pygame.display.set_mode((Width, Height), 0, 32)
# #    screen = display.set_mode([Width, Height])
#     pygame.display.set_caption("Activity")
#     bg = 0
#     screen.fill(bg)
#     # while True:
#     #     pygame.draw.rect(screen, [10, 10, 10],(0, 0 , Width, Height + 1), 0)
#     #     pygame.display.update()

    argv = sys.argv
    argv = argv[1:]
    sys.stderr.writelines(str(argv)+'\n')
    proc = sp.Popen(argv,
                    shell=True,
                    stdout = sp.PIPE)
    while True:
        dat = proc.stdout.readline()
        print dat,


def main():
    global ymax

    argv = sys.argv
    argv = argv[1:]
    sys.stderr.writelines(str(argv)+'\n')
    proc = sp.Popen(argv,
                    stdout = sp.PIPE)
    

    # foo = Graph(0,250,Width,Height,"Activity")
    screen = display.set_mode([Width, Height])
    display.set_caption("Activity")
    draw.rect(screen, [10, 10, 10],(0, 0 , Width, Height + 1), 0)
    display.update()
    actcnt = 0;
    colors = {}
    cnt = 0;
    actpre = None
    while True:
 #       sys.stderr.write('---\n')
        # get input
        for ee in event.get():
            if ee.type == QUIT:
                return
            elif ee.type == KEYDOWN:
                if ee.key == K_ESCAPE:
                    return
                elif ee.key == K_PLUS:
                    ymax = ymax * 2
                    msg = 'new ymax ='+str(ymax)+'\n'
                    sys.stderr.write(msg)
                elif ee.key == K_KP_PLUS:
                    ymax = ymax * 2
                    print 'new ymax =',ymax
                elif ee.key == K_EQUALS:
                    ymax = ymax * 2
                    print 'new ymax =',ymax
                elif ee.key == K_MINUS:
                    ymax = ymax / 2
                    print 'new ymax =',ymax
            """
            yval = cnt%Height
            col = palette[cnt%255]
            trace(screen,[(col,yval)])
            cnt += 1
            """
#        try:
        dat = proc.stdout.readline()
#        sys.stderr.writelines('---|||'+dat)
#        except:
#            print 'had a problem reading activity pipe.'
        dat = dat.split()
# dat has (string,count) pairs on each line.
#        dat = [int(x) for x in dat.split()]
        ldat = len(dat)
        if ldat == 0:
            return
        if ldat % 2:
            print "length of activity output not multiple of 2!"
            print '++++++++++ '
            print dat
            print '---------------\n'
        if datfile != None:            # write to external file
            for dd in dat:
                datout.write(str(dd)+' ')
            datout.write("\n")
        activity = {}
        actflat = {}
        fofo = range(0,ldat,2)          # swizzle data to graph: idx of beginning of each data chunk
        for i in range(len(fofo)):
            activity[dat[fofo[i]]] = int(dat[fofo[i]+1])
        dat = [(dat[fofo[i]],int(dat[fofo[i]+1])) for i in range(len(fofo))]
        # dat = [(x,y) for [x,y] in dat]

        for xx in activity:
            if xx not in colors:
                if cnt==0:              # 1st data chunk all white
                    colors[xx] = (255,255,255)
                else:
                    colors[xx] = palette[actcnt%255]
                    actcnt += 1
        if actpre is not None:
            for xx in activity:
                if xx in actmax:
                    if activity[xx] <= actmax[xx]:   # flat detected
                        activity[xx] = None
                    else:
                        actmax[xx] = activity[xx]
        trdat = [(colors[xx],activity[xx]) for xx in activity if activity[xx] is not None]   # assemble data for trace
        if len(trdat) > maxact:
            trdat = trdat[0:maxact]
        trace(screen,trdat)
        # if datfile != None:            # write to external file
        #     for xx,act in dat:
        #         datout.write(str(act)+' ')
        #     datout.write("\n")
        cnt += 1                        # counting number of data chunks piped from source program

if __name__=='__main__':
    main()

