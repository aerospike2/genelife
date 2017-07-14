#!/usr/bin/env python

import pygame, sys
from pygame.locals import *

import random

from numpy.random import randint
from numpy.random import rand
import numpy as np

import time

class VideoData:
    """ sample video data source"""

    #def getpixels(self,ix,iy):
    #    return self.rnr*256*256*256 + 2*iy*256 + 2*ix*256*256
    
    def getframe(self):
        self.rnr = randint(0,0xff)
        for ix in range(video.W):
            for iy in range(video.H):
                video.pix[ix][iy] = self.rnr*256*256*256 + 2*iy*256 + 2*ix*256*256

    def __init__(self):
        self.frame_function = self.getframe
        #self.pixel_function = self.getpixels

class FastVideo:
    """Fast video graphics for 128*128 array using pygame"""

    def __init__(self,xself=None,dt=0.0,nframes=0,Width=128,Height=128):
        self.dt=dt
        self.xself=xself
        self.nframes=nframes
        self.Seed = 0       # random seed
        pygame.init()
        self.step=0
        random.seed(self.Seed)
        self.W = Width
        self.H = Height
        # print self.W,self.H
        self.surf = pygame.Surface((Width,Height),0,32) # working surface
        self.pix = pygame.PixelArray(self.surf)     # pixel array
        self.surf2 = pygame.Surface((2*self.W,2*self.H),0,32)       # intermediate surface for 4x zoom
        self.surfdisp = pygame.display.set_mode((4*self.W, 4*self.H), 0, 32) #  this defines the surface (zoomed) for the display window
        bg = 0                          # color black
        self.surf.fill(bg)
        self.disp()

    def disp(self):
        # set up 2x and 4x surfaces for 2x * 2x zoom just before display update
        pygame.transform.scale2x(self.surf,self.surf2)
        pygame.transform.scale2x(self.surf2,self.surfdisp)
        pygame.display.update()
    
    def update(self,xself):
        xself.frame_function()
        #for ix in range(self.W):
        #    for iy in range(self.H):
        #        self.pix[ix][iy] = xself.pixel_function(ix,iy)
        self.disp()   # 2x * 2x zoom

        self.step = self.step+1
        if self.step == self.nframes:
            return False
        for event in pygame.event.get():
            if event.type == QUIT:
                return False
        return True


    def run(self,xself):
        xself.play = True
        t0 = time.time()
        while xself.play is True:
            xself.play = self.update(xself)
            t1 = time.time()
            t1m0=t1-t0
            if (self.dt - (t1m0)) < 0.0 :
                error = "overtime"
            else:
                error = ""
                time.sleep(self.dt-(t1m0))
            t0 = time.time()
            print self.step,time.ctime(),'delta time = ',t1m0,error,self.dt
            #xself.top.update_idletasks()
        pygame.quit()
 
 
if __name__ == "__main__":
    video_src = VideoData()
    video = FastVideo(video_src,0.017,100)
    video.run(video_src)
