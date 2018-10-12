#!/usr/bin/env python

import pygame, sys
from pygame.locals import *
import random
from numpy.random import randint
from numpy.random import rand
import time
import string
import re

##################################################
## configuration params

# world
Height = 64
Width = 64
Dens = 0.01                                # initial density of agents

# agent ID specification
LenID = 64    # length of ID string
Nletters = 8   # number of letters used  e.g. 4 => letters drawn from ABCDabcd
Nrelevant = 2       # number of relevant letters (others are neutral)
Seed = 0       # random seed

# reproduction
Mutrate = 0.02
Duprate = 0.01
Dierate = 0.001

## configuration params
##################################################


random.seed(Seed)

def ranstr(size=6, chars=string.ascii_uppercase + string.ascii_lowercase):
    return ''.join(random.choice(chars) for _ in range(size))


def ranchar(chars = string.ascii_uppercase + string.ascii_lowercase):
    return random.choice(chars)

# setup letter universe
Letters = ''.join(string.ascii_uppercase[i]+string.ascii_lowercase[i] for i in range(Nletters))
Relevant = ''.join(string.ascii_uppercase[i]+string.ascii_lowercase[i] for i in range(Nrelevant))
Irrelevant = ''.join([x for x in Letters if x not in Relevant])
Relattack = ''.join(string.ascii_uppercase[i] for i in range(Nrelevant))
Reldefend = ''.join(string.ascii_lowercase[i] for i in range(Nrelevant))


def findattacks(str):                      # find and return attack strings
    re_attack = '['+Relattack+']*'
    ret = re.findall(re_attack,str)
    ret = [ret[i] for i in range(len(ret)) if ret[i] is not '']
    return ret

def finddefends(str):                      # find and return attack strings
    re_attack = '['+Reldefend+']*'
    ret = re.findall(re_attack,str)
    ret = [ret[i] for i in range(len(ret)) if ret[i] is not '']
    return ret

#  choose random chunk of parent and duplicate in random place of child
def dupgene(id1,id2):                     # id1 from parent, id2 from child
# commented out prints are for validating and testing this function...
    assert len(id1) == len(id2)
    id = id1
    tmp1 = random.choice(range(len(id)))
    tmp2 = random.choice(range(len(id)))
    start = min([tmp1,tmp2])
    end = max([tmp1,tmp2])
    gene = id[start:end]
#    print 'start, end = ',start,end
#    print 'gene = ',gene
    targ = random.choice(range(len(id2)))
#    print 'target = ',targ
    if targ + len(gene) > len(id):    # need to trim...
        cut = targ+len(gene) - len(id)
        gene = gene[0:(len(gene)-cut)]
#        print 'new gene = ',gene
    ret = id2[:targ] + gene + id2[(targ+len(gene)):]
    return ret



class Agent:                              #  agent will have (x,y,ID_vector)
    def __init__(self,x,y,id):
        self.x = x
        self.y = y
        self.id = id
        refix = '['+Irrelevant+']'
        self.idrel = re.sub(refix,'-',id) # id collapsed to only relevant letters
        self.popidx = None                # index in population agent array
        self.attacks = findattacks(id)
        self.defends = finddefends(id)
        self.color = self.Color()         # hex string to color pixels
        self.activity = 0

    def Color(self,scale=40):
        if len(self.attacks) > 0: 
            Natt = max([len(x) for x in self.attacks])
        else:
            Natt = 0
        if len(self.defends) > 0:
            Ndef = max([len(x) for x in self.defends])
        else:
            Ndef = 0
        Natt = min(scale*Natt,255)
        Ndef = min(scale*Ndef,255)
        col = '0x'+''.join('%02x'%i for i in [Natt,Ndef,0]) # form the hex string for the color
        return col

class Latmap:

    def __init__(self, dens=0.5):
        pygame.init()
        self.surf = pygame.Surface((Width,Height),0,32)
        self.pix = pygame.PixelArray(self.surf)     # pix[x][y] = color 
        self.lattice = [[None for x in range(Width)] for y in range(Height)]
        self.dens = dens
        self.agents = []
        self.init_agents()                # set of randomly placed Agents w/ random ID

    def countlattice(self):
        cnt = 0
        for x in range(Width):
            for y in range(Height):
                if self.lattice[x][y] is not None:
                    cnt += 1
        return cnt

    def consistency(self):
        for a in self.agents:
            if a is not None:
                if self.lattice[a.x][a.y] != a.popidx:
                    print 'prob at ',a.x,a.y,': popidx ',a.popidx
                    return False
        return True

    def nodup(self):
        idx = []
        for a in self.agents:
            if a is not None:
                idx.append(a.popidx)
        if len(idx) != len(set(idx)):
            return False
        else:
            return True

    def install_agent(self,aa):
        aa.popidx = len(self.agents)      # index of this new agent in agents list
        self.agents.append(aa)
        self.lattice[aa.x][aa.y] = aa.popidx
        self.pix[aa.x][aa.y] = pygame.Color(aa.color)
# for debugging...
#        assert self.nodup()
#        assert self.consistency()
            

    def uninstall_agent(self,aa):
        self.agents[aa.popidx] = None
        self.lattice[aa.x][aa.y] = None
        self.pix[aa.x][aa.y] = 0        # black background...

    def agent_exists(self,aa):            #for debugging...
        if self.agents[aa.popidx] == None:
            return False
        if self.lattice[aa.x][aa.y] == None:
            return False
        return True
        

    def compactify_agents(self):          # compactify agent list after some deaths
#        print 'Before compactify  Nagents vs. countlattice:  ',len(self.agents),self.countlattice()
        self.agents = [self.agents[i] for i in range(len(self.agents)) if self.agents[i] is not None]
        for i in range(len(self.agents)):   # reset popidx for all agents...
            self.agents[i].popidx = i
            x = self.agents[i].x
            y = self.agents[i].y
            self.lattice[x][y] = i
#        print 'After compactify Nagents vs. countlattice:  ',len(self.agents),self.countlattice()

    def init_agents(self):
        N = int(self.dens*Width*Height)
        allxy = [(x,y) for x in range(Width) for y in range(Height)]
        random.shuffle(allxy)             # randomize list of xy pairs
        # print 'initializing',N,'agents...'
        for i in range(N):
            x, y = allxy[i]               # a random xy pair
            id = ranstr(LenID,Letters)
            aa = Agent(x,y,id)
            self.install_agent(aa)
        # print 'done.'
        
    # predator action:
    # go through all attack strings, check to see if prey defends each one
    # if yes for any, no action.
    # if no, prey is disappears and predator reproduces.
    def predate(self,a,aa):
        assert self.agent_exists(a)
        assert self.agent_exists(aa)
        assert a.popidx != aa.popidx      # must be different agents
        attacks = a.attacks
        defends = aa.defends
        # for neutral model, choose attack result randomly:
        # if rand() > 0.5:
        #     attackwon = 1
        # else:
        #     attackwon = 0
        attackwon = 0                     # default for no attacks
        for x in attacks:
            attackwon = 1
            for y in defends:
                yy = y.upper()
                if x in yy:
                    attackwon = 0         # successful defense
                    break;                # defense succeeded, go to next attacker.
        if attackwon:                     # no successful defense
            x = aa.x                      # grab coordinates of prey
            y = aa.y
            self.uninstall_agent(aa)      # remove unsuccessful prey
            self.reproduce(a,x,y)              # reproduce attacker in place of prey
                    
    def reproduce(self,a,x,y):
        assert self.lattice[x][y] == None
        assert self.agent_exists(a)
        id = a.id
        ranums = rand(LenID)
        for i in range(len(ranums)):
            if ranums[i] < Mutrate:
                id = id[:i] + ranchar(Letters) + id[(i+1):]
        aa = Agent(x,y,id)           # also recomputes color, attackers and defenders
        ## now do 'gene duplication'...
        if rand() < Duprate:
            aa.id = dupgene(a.id,aa.id)
        # end duplication...
        self.install_agent(aa)           # 
        
    def update(self):
        w = Width
        h = Height
        popdirty = 0                      # to keep track of removals to compact agent list at end
        idx = range(len(self.agents))
        random.shuffle(idx)               # for random updating of all agents
        for ii in idx:
            a = self.agents[ii]
            if a is None:                 # oops! hit one that died earlier this update
                continue
            self.agents[ii].activity += 1
            x = a.x
            y = a.y
            if rand()<Dierate:
                self.uninstall_agent(a)
                continue
            l = (x - 1 + w) % w
            r = (x + 1) % w
            u = (y + 1 ) % h
            d = (y - 1 +h) % h 
            nbrs = ((l,y),(r,y),(x,u),(x,d))
            # check 2 random nbrs...
            nbridx = range(len(nbrs))
            random.shuffle(nbridx)
            for nbr in [nbrs[nbridx[0]],nbrs[nbridx[1]]]:
                xx,yy = nbr
                aidx = self.lattice[xx][yy]
                if aidx is not None:
                    anbr = self.agents[aidx]
                    assert self.agent_exists(anbr)
                    assert self.agent_exists(a)
                    if a.popidx == anbr.popidx:   # error: duplicate
                        # print 'a xy, popidx = ',a.x,a.y,a.popidx
                        # print 'anbr xy, popidx = ',anbr.x,anbr.y,anbr.popidx
                        print nbrs
                        print 'cur nbr = ',nbr
                        assert False
                    self.predate(a,anbr)
                else:
                    self.reproduce(a,xx,yy)
        # painful but got to compactify the list when one is removed:
        self.compactify_agents()
# for QA... leave out for speed.
#        assert self.consistency()
#        assert self.nodup()
            
    def run(self, sleep=0):
        w, h = Width, Height
        agents = self.agents
 
        
        surfdisp = pygame.display.set_mode((4*w, 4*h), 0, 32)
        surf2 = pygame.Surface((2*Width,2*Height),0,32)
        # surf set in __init__
        bg = 0                          # color black
        self.surf.fill(bg)
 
        step = 0
        event = None
        pygame.transform.scale2x(self.surf,surf2)
        pygame.transform.scale2x(surf2,surfdisp)
        pygame.display.update()
        # print 'entering event loop'
#        fact = open('/tmp/activity','w')
        idact = {}                        # dictionary to hold activities
        while event is None or event.type != QUIT:
            # sys.stderr.writelines(str(step)+'   '+str(len(self.agents))+'\n')
            step += 1
            self.update()
            ####################################################
            ## here is the activity computation 
            for a in self.agents:
                for x in a.attacks:
                    if x in idact:
                        idact[x] += 1
                else:
                    idact[x] = 1
                for x in a.defends:
                    if x in idact:
                        idact[x] += 1
                else:
                    idact[x] = 1
            strout = ''
            for x in idact:
                strout = strout + x + ' ' + str(idact[x]) + ' '
            print strout
            sys.stdout.flush()
            ## end activity computation 
            ####################################################
            pygame.transform.scale2x(self.surf,surf2)
            pygame.transform.scale2x(surf2,surfdisp)
            pygame.display.update()
            for event in pygame.event.get():
                if event.type == QUIT:
                    break
                if sleep > 0:
                    time.sleep(sleep)
    pygame.quit()
 
 
if __name__ == "__main__":
    latmap = Latmap(dens=Dens)
    latmap.run()
