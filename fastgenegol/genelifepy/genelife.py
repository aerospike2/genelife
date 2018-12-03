import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpcolors
from   matplotlib.colors import ListedColormap
import matplotlib.animation as animation
import matplotlib
import genelife_update_module as genelife
import pygame as pg


log2N = genelife.get_log2N()                       #  log2N=7 => N=128 get value from C library which is where it must be changed
N = 2**log2N
N2 = N*N
Nmask = N-1
NbP = 1
gol = np.zeros(N2,np.uint64)
golg = np.zeros(N2,np.uint64)
golgstats = np.zeros(N2,np.uint64)
                                    # graphics
cgrid = np.zeros((N,N),np.int32)
cgolg =np.zeros(N2,np.int32)
colorfunction = 0
scr= None
screen = None
scalex2 = False
cancol=[]
caption = None
selectiontext0007 = ["largest value","most ones","scissors-well-stone-paper","not well ordered","two target","predator prey","cooperative","neutral"];
selectiontext0815 = ["sum fixed","sum variable","edge fixed","edge variable","canonical fixed","canonical variable","match 2nd layer","NYI"];
selectiontext1623 = ["2-16 plane pairwise","2-16 plane pairwise","2-16 plane nearby","2-16 plane nearby","2-64 plane matching","2-64 plane matching","2-64 plane matching","2-64 plane matching"]

Width = N                           # value specified in imported genelife.py
Height = N
dispinit = False
updatesenabled = True
displayplanes=0xffff
displayoneplane=64
mat = []
ndisp = 100
nskip = 0
nhist = 0
nstat = 0

gogo = True
pixeldat = ""
paramdat = ""
mouseclicked = False
mouseclicked2 = False
pause = 0
ymax = 10000
maxPlane = 4
offdx = offdy=offdt=0
quadrants = -1
oldymax = ymax

                                    # counter and toggle initialization
nrun = 1
cnt = 0
framenr = 0
savecnt = 0                         # counter for saved images
randomsoup = 0
vscrolling = 0
                                    # parameter initialization
runparams = np.zeros(7,np.int32)    # 7 parameters passed to C
simparams = np.zeros(5,np.int32)    # 5 parameters passed to C

selection = 0
rulemod = 0
repscheme=0
survivalmask=0
overwritemask=0
surviveover = np.array([survivalmask,overwritemask],dtype=np.uint32)

ncoding=0
nlog2pmut=0
initfield=100
startgenechoice=0
initial1density=16284
initialrdensity = 32768
                                    # offset initialization
offsets = [[-1, 0,-1],
           [ 1, 0,-1],
           [ 0,-1,-1],
           [ 0, 1,-1]]
numHis = len(offsets)
histo=np.zeros(numHis,np.uint64)
flatoff =  [x for sublist in offsets for x in sublist]
npoffsets = np.array(flatoff,np.int32)

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

    if type not in ('bright', 'soft', 'grad'):
        print ('Please choose "bright" or "soft" or "grad" for type')
        return
    if verbose:
        print('Number of labels: ' + str(nlabels))
    np.random.seed(123456)
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
#            randRGBcolors[1] = [0, 0, 1] # this color otherwise is black
        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    if type == 'grad':
#        randRGBcolors = np.array((np.linspace(0.0, 1.0, num=256),
#                                  np.linspace(1.0, 0.0, num=256),
#                                  np.concatenate(np.linspace(1.0, 0.0, num=128),np.linspace(0.0, 1.0, num=128)))).T
        randRGBcolors =[((i-1)/255.,(255.-(i-1))/255.,(i-1)/128. if i<=128 else (255.-(i-1))/128.) for i in xrange(nlabels)]
        randRGBcolors[0] = [0, 0, 0]
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

def colorgrid():
    """ colors array according to grid and genegrid using colormethod"""
    global cgrid,cgolg,crid2,N
    genelife.colorgenes(cgolg)
    cgrid2=np.reshape(cgolg,(N,N)).T
    cgrid[:,0:N] = cgrid2   # is there a faster version of this copy that moves the data?
    return
#-----------------------------------------------------------------------------------------------------------

def packrepscheme(repscheme,survivalmask,overwritemask):
    if survivalmask<4 and overwritemask<4:
        repscheme = repscheme + (survivalmask<<24) + (overwritemask<<26)
    else:
        print "Error: can't pack masks, they are too large!"
    return(repscheme)
#-----------------------------------------------------------------------------------------------------------

def init_button_arrays():
    """ initialize information for button area at base of display"""
    global ncanon,cancolors,cancol
                                                      # initialize lists to empty
    ncanon = []                                       # lists of number of successive buttons with the same color
    cancolors = []                                    # lists of button colors
    cancol =[]                                        # expanded lists of button colors, by button
                                                      # no of entries per color region
    ncanon.append([2,2,2,2,2,4,2,2])                  # selection 0-7
    ncanon.append([8,8])                              # selection 8,9
    ncanon.append([3,4,5,4,3,3,4,5,4,3])              # selection 10,11
    ncanon.append([4,7,10,7,4,4,7,10,7,4])            # selection 12,13
    ncanon.append([16,1,1,1,1,1])                     # selection 16-19
                                                      # colors for different color regions for buttons, colorvals must be < 128
    cancolors.append([[0,100,0],[0,50,100],[0,80,80],[100,0,0],[100,100,0],[50,100,0],[0,0,127],[0,100,50]]) # selection 0-7
    cancolors.append([[0,0,127],[0,100,0]])          # selection 8,9
    cancolors.append([[100,0,0],[0,100,0],[0,0,127],[0,100,0],[100,0,0],[100,0,0],[0,100,0],[0,0,127],[0,100,0],[100,0,0]]) # selection 10,11
    cancolors.append([[100,0,0],[0,100,0],[0,0,127],[0,100,0],[100,0,0],[100,0,0],[0,100,0],[0,0,127],[0,100,0],[100,0,0]]) # selection 12,13
    cancolors.append([[100,0,0],[100,100,0],[0,80,80],[0,0,127],[0,100,0],[80,0,80]]) # selection 16-19
                                                     # lists of colors for individual buttons expanded from above, first initialize to zero
    cancol.append(np.zeros((18,3),np.int32))
    cancol.append(np.zeros((16,3),np.int32))
    cancol.append(np.zeros((38,3),np.int32))
    cancol.append(np.zeros((64,3),np.int32))
    cancol.append(np.zeros((21,3),np.int32))
    
    for l in range(len(ncanon)):                    # buttons for different selection schemes
        k=0
        for j in range(len(ncanon[l])):
            for i in range(ncanon[l][j]):
                cancol[l][k]=cancolors[l][j]
                k = k+1
#-----------------------------------------------------------------------------------------------------------
def init_buttons():    # initialize parameter buttons
    global repscheme,survivalmask,overwritemask,selection,ncoding,displayplanes
    global cancol
    global scr
    global Height,Width
    global log2N,NbP

    init_button_arrays()
    pg.draw.rect(scr,[50,50,50],[0,Height+8,Width,5])
    if selection<8:
        for k in range(18):
            if k<14:
                bit = (repscheme>>k)&0x1
            elif k<16:
                bit = (survivalmask>>(k-14))&0x1
            elif k<18:
                bit = (overwritemask>>(k-16))&0x1
            pg.draw.rect(scr,cancol[0][k]*(1+bit),[k<<(log2N-6),Height+8,3,5])
    elif selection<10:
        for k in range(16):
            pg.draw.rect(scr,cancol[1][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+8,3,5])
    elif selection<12:
        for k in range(19):
            pg.draw.rect(scr,cancol[2][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+8,3,5])
            pg.draw.rect(scr,cancol[2][k+19]*(1+((overwritemask>>(k))&0x1)),[(k+19)<<(log2N-6),Height+8,3,5])  
    elif selection<14:
        for k in range(32):
            pg.draw.rect(scr,cancol[3][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+8,3,5])
            pg.draw.rect(scr,cancol[3][k+32]*(1+((overwritemask>>(k))&0x1)),[(k+32)<<(log2N-6),Height+8,3,5])               
    elif selection>=16 and selection<=19:
        NbP = (ncoding>>16)&0xf
        displayplanes=(0x1<<NbP)-1
        if not NbP:
            NbP = 16
        for k in range(21):
            if k<NbP:
                pg.draw.rect(scr,cancol[4][k]*2,[k<<(log2N-6),Height+8,3,5])
            elif k<16:
                pg.draw.rect(scr,[80,80,80],[k<<(log2N-6),Height+8,3,5]) // grey
            elif k<21:
                bit = (repscheme>>(k-16))&0x1
                pg.draw.rect(scr,cancol[4][k]*(1+bit),[k<<(log2N-6),Height+8,3,5])
#-----------------------------------------------------------------------------------------------------------

def display_init():
    global screen,scr,scalex2,Width,Height
    global caption,cnt,cgrid,dispinit

    dispinit = True
    if (Height <=512):
        scalex2 = True
        screen = pg.display.set_mode([2*Width, 2*Height+32])     # opens the pygame window
        scr = pg.surface.Surface([Width,Height+16], 0)
    else:
        scalex2 = False
        screen = pg.display.set_mode([Width, Height])

    cnt = 0
    caption = "Gene Life at iteration %d" % cnt
    pg.display.set_caption(caption)
    
    if scalex2:
        pg.transform.scale2x(scr,screen)
        pg.display.update()
        cgrid = pg.surfarray.pixels2d(scr)
    else:
        pg.display.update()
        cgrid = pg.surfarray.pixels2d(screen)
#-----------------------------------------------------------------------------------------------------------

def show0(count=True):
# display initial population and count species
    global framenr
    global scr, screen, scalex2,caption
    global repscheme,survivalmask,overwritemask,selection
    global cancol
    global dispinit
    
    if not dispinit:
        display_init()
    caption = "Gene Life at iteration %d" % framenr
    pg.display.set_caption(caption)

    init_buttons()                           # initialize parameter buttons
    
    colorgrid()
    # pg.transform.scale2x(scr,screen)       # use this for standard dithered display
    if scalex2:
        pg.transform.scale2xact(scr,screen)  # use this for custom pygame no smoother such as in scale2x
    pg.display.flip()
    if(count):
        genelife.countspecieshash()
#-----------------------------------------------------------------------------------------------------------

def step(count=True):
    """single step and update display and species counts"""
    global framenr
    #global gol,golg,golgstats
    global scr, screen, scalex2
    global dispinit
    
    if not dispinit:
        display_init()

    update_sim(count)
    caption = "Gene Life at iteration %d" % framenr
    pg.display.set_caption(caption)

    # pg.transform.scale2x(scr,screen)  # use this for standard dithered display
    if scalex2:
         pg.transform.scale2xact(scr,screen)  # use this for custom pygame no smoother
    pg.display.flip()
    if (count):
        genelife.countspecieshash()
#-----------------------------------------------------------------------------------------------------------

# infinite loop of display updates
# cmd click in graphics window to stop, click for pixel details or quadrant selection,
# alt-click or arrow keys for recolor,
# +/- keys reserved for activity ymax : actually the crossover value in N* act/(ymax+act)
# keys lower case - decrement, upper case - increment, alt - input value: y,Y ymax q,Q quadrant
# misc. keys save image
def run(count=True):
    global framenr
    global scr, screen, scalex2
    global N,NbP
    global gol,golg,golgstats
    global colorfunction
    global ymax
    global updatesenabled
    global rulemod,repscheme,survivalmask,overwritemask,selection,ncoding,displayplanes
    global savecnt
    global cancol
    global Height,Width
    global dispinit
    global randomsoup,vscrolling
    global selectiontext0007, selectiontext0815, selectiontext1623
    global surviveover
    global gogo,pause,mouseclicked,mouseclicked2,pixeldat,paramdat
    global ymax,maxPlane,offdx,offdy,offdt,quadrants,oldymax,displayoneplane
    global parhelp
    
    if not dispinit:
        display_init()
    init_buttons()
    
    surviveover = np.array([survivalmask,overwritemask],dtype=np.uint32)
    gogo = True
    pixeldat = ""
    paramdat = ""
    mouseclicked = False
    mouseclicked2 = False
    pause = 0
    ymax = 10000
    maxPlane = 4
    offdx = offdy=offdt=0
    quadrants = -1
    oldymax = genelife.setget_act_ymax(ymax)
    displayoneplane=64

    if selection>=16 & selection<19:
        displayplanes = (0x1<<NbP)-1

    while (gogo):
        for event in pg.event.get():
            if event.type==pg.QUIT:
                mouseclicked = False
                gogo = False
            if event.type==pg.MOUSEBUTTONDOWN:
                if event.button == 2:          # quit event loop on middle mouse button (option-click)
                    mouseclicked = False
                    gogo = False
                elif event.button == 1:          # get mouse coords on mouse event
                    mouseclicked = True
                    mouse_pos = pg.mouse.get_pos()
                    x = (int) (mouse_pos[0]//2)
                    y = (int) (mouse_pos[1]//2)
                    if y >= N:
                        k=x>>(log2N-6)
                        if selection<8:
                            if k<18:
                                if k<14:
                                    repscheme = repscheme ^ (1<<k)
                                    bit = (repscheme>>k)&0x1
                                    print ("repscheme changed to %x" % (repscheme))
                                elif k<16:
                                    survivalmask = survivalmask ^ (1<<(k-14))
                                    bit = (survivalmask>>(k-14))&0x1
                                    print ("survivalmask changed to %x" % (survivalmask))
                                else:
                                    overwritemask = overwritemask ^ (1<<(k-16))
                                    bit = (overwritemask>>(k-16))&0x1
                                    print ("overwritemask changed to %x" % (overwritemask))
                                survivalmask
                                pg.draw.rect(scr,cancol[0][k]*(1+bit),[k<<(log2N-6),Height+8,3,5])
                                surviveover[0],surviveover[1]= survivalmask,overwritemask
                                genelife.set_surviveover64(surviveover)
                                genelife.set_repscheme(repscheme)
                        elif selection < 10:
                            if k<16:
                                survivalmask = survivalmask ^ (1<<k)
                                print ("survivalmask changed to %x" % (survivalmask))
                                pg.draw.rect(scr,cancol[1][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+8,3,5])
                                surviveover[0],surviveover[1]= survivalmask,overwritemask
                                genelife.set_surviveover64(surviveover)
                        elif selection < 12:
                            if k<38:
                                if k<19:
                                    survivalmask = survivalmask ^ (1<<k)
                                    print ("survivalmask changed to %x" % (survivalmask))
                                    pg.draw.rect(scr,cancol[2][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+8,3,5])
                                else:
                                    overwritemask = overwritemask ^ (1<<(k-19))
                                    print ("overwritemask changed to %x" % (overwritemask))
                                    pg.draw.rect(scr,cancol[2][k]*(1+((overwritemask>>(k-19))&0x1)),[k<<(log2N-6),Height+8,3,5])
                                surviveover[0],surviveover[1]= survivalmask,overwritemask
                                genelife.set_surviveover64(surviveover)
                        elif selection<14:
                            if k<64:
                                if k<32:
                                    survivalmask = survivalmask ^ (1<<k)
                                    print ("survivalmask changed to %x" % (survivalmask))
                                    pg.draw.rect(scr,cancol[2][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+8,3,5])
                                else:
                                    overwritemask = overwritemask ^ (1<<(k-32))
                                    print ("overwritemask changed to %x" % (overwritemask))
                                    pg.draw.rect(scr,cancol[2][k]*(1+((overwritemask>>(k-32))&0x1)),[k<<(log2N-6),Height+8,3,5])
                                surviveover[0],surviveover[1]= survivalmask,overwritemask
                                genelife.set_surviveover64(surviveover)
                        elif selection < 16:
                                pass
                        elif selection < 20:
                            if k<NbP:
                                displayplanes = displayplanes ^ (1<<k)
                                pg.draw.rect(scr,cancol[3][k]*(1+((displayplanes>>k)&0x1)),[k<<(log2N-6),Height+8,3,5])
                                genelife.set_displayplanes(displayplanes)
                            elif k>=16 and k<21:
                                repscheme = repscheme ^ (1<<(k-16))
                                print ("repscheme changed to %x" % (repscheme))
                                pg.draw.rect(scr,cancol[3][k]*(1+((repscheme>>(k-16))&0x1)),[k<<(log2N-6),Height+8,3,5])
                                genelife.set_repscheme(repscheme)
                    else: # selection >= 20
                        if colorfunction < 4 or colorfunction == 8:
                            genelife.get_curgol(gol)    # get current gol,golg,golgstats arrays
                            genelife.get_curgolg(golg)
                            genelife.get_curgolgstats(golgstats)
                            if quadrants >= 0:   # set the two bits in repscheme corresponding to quadrant
                                repscheme=genelife.set_repscheme_bits(quadrants,x,y,surviveover)
                                survivalmask  = surviveover[0]
                                overwritemask = surviveover[1]
                                repscheme=packrepscheme(repscheme,survivalmask,overwritemask)
                                print ("repscheme changed to %x" % (repscheme))
                                quadrants = -1
                                pixeldat = ""
                            else:
                                pixeldat = "(%d,%d) gol %016x gene %016x status %016x" % (x,y,gol[x+y*N],golg[x+y*N],golgstats[x+y*N])
                        elif colorfunction == 4:
                            genelife.get_acttrace(golg)
                            pixeldat = "(%d,%d) gene %016x" % (x,y,golg[x+y*N])
                        elif colorfunction <= 7:
                            genelife.get_genealogytrace(golg)
                            pixeldat = "(%d,%d) gene %016x" % (x,y,golg[x+y*N])
                            genelife.set_selectedgene(golg[x+y*N])
                            print framenr,pixeldat
                elif event.button == 3:          # single plane choice right mouse button (-click)
                    mouse_pos = pg.mouse.get_pos()
                    x = (int) (mouse_pos[0]//2)
                    y = (int) (mouse_pos[1]//2)
                    if y >= N and selection>=20:
                        k=x>>(log2N-6)
                        if k<64:
                            displayoneplane=k
                            genelife.set_displayoneplane(displayoneplane)
                    mouseclicked2 = True

            elif event.type==pg.MOUSEBUTTONUP:
                mouseclicked = False
                mouseclicked2 = False
                if selection>=20:
                    displayoneplane=64
                    genelife.set_displayoneplane(displayoneplane)
                    if not updatesenabled:
                        updatesenabled=True
            elif event.type==pg.MOUSEMOTION:
                if mouseclicked:
                    mouse_pos = pg.mouse.get_pos()
                    x = mouse_pos[0]//2
                    y = mouse_pos[1]//2
                    if x < N and y < N:
                        if colorfunction < 4 or colorfunction == 8:
                            genelife.get_curgol(gol)    # get current gol,golg,golgstats arrays
                            genelife.get_curgolg(golg)
                            genelife.get_curgolgstats(golgstats)
                            pixeldat = "(%d,%d) gol %016x gene %016x status %016x" % (x,y,gol[x+y*N],golg[x+y*N],golgstats[x+y*N])
                        elif colorfunction == 4:
                            genelife.get_acttrace(golg)
                            pixeldat = "(%d,%d) gene %016x" % (x,y,golg[x+y*N])
                        elif colorfunction <= 7:
                            genelife.get_genealogytrace(golg)
                            pixeldat = "(%d,%d) gene %016x" % (x,y,golg[x+y*N])
                            genelife.set_selectedgene(golg[x+y*N])
                elif mouseclicked2:
                    if colorfunction == 2:
                        mouse_pos = pg.mouse.get_pos()
                        x = (int) (mouse_pos[0]//2)
                        y = (int) (mouse_pos[1]//2)
                        if y >= N and selection==12:
                            k=x>>(log2N-6)
                            if k<64:
                                displayoneplane=k
                                genelife.set_displayoneplane(displayoneplane)
                            if updatesenabled:
                                updatesenabled=False
            elif event.type == pg.KEYDOWN:
                if event.key == pg.K_h:
                    parhelp()
                elif event.key == pg.K_SPACE:
                    pause = (pause+1)%2
                elif event.key == pg.K_RIGHT:
                    colorfunction = (colorfunction + 1) % 9
                    genelife.set_colorfunction(colorfunction)
                elif event.key == pg.K_LEFT:
                    colorfunction = (colorfunction - 1) % 9
                    genelife.set_colorfunction(colorfunction)
                elif event.key == pg.K_PLUS or event.key == pg.K_KP_PLUS or event.key == pg.K_EQUALS:
                    ymax = ymax * 2
                    oldymax = genelife.setget_act_ymax(ymax)
                    print 'new ymax =',ymax
                elif event.key == pg.K_MINUS:
                    ymax = ymax / 2
                    oldymax = genelife.setget_act_ymax(ymax)
                    print 'new ymax =',ymax
                elif event.key == pg.K_x:
                    if pg.key.get_mods() & pg.KMOD_SHIFT: offdx = offdx+1
                    else: offdx = offdx-1
                    print "offset dx changed to ",offdx
                    genelife.set_offsets(offdx,offdx,offdt)
                elif event.key == pg.K_y:
                    if pg.key.get_mods() & pg.KMOD_SHIFT: offdy = offdy+1
                    else: offdy = offdy-1
                    print "offset dy changed to ",offdy
                    genelife.set_offsets(offdx,offdy,offdt)
                elif event.key == pg.K_t:
                    if pg.key.get_mods() & pg.KMOD_SHIFT:
                        if(offdt<0): offdt = offdt+1
                    elif offdt>-maxPlane+1: offdt = offdt-1
                    print "offset dt changed to ",offdt
                    genelife.set_offsets(offdx,offdy,offdt)
                elif event.key == pg.K_q:
                    if pg.key.get_mods() & pg.KMOD_ALT:
                        quadrants = input("Enter an integer between -1 and 6: ")
                    elif pg.key.get_mods() & pg.KMOD_SHIFT:
                        if quadrants < 7: quadrants = quadrants+1
                    else:
                        if quadrants > 0: quadrants = quadrants-1
                    print "quadrants changed to ",quadrants
                    genelife.set_quadrant(quadrants)
                elif event.key == pg.K_r:
                    if pg.key.get_mods() & pg.KMOD_SHIFT:
                        randomsoup=1-randomsoup
                        print "randomsoup changed to ",randomsoup
                        genelife.set_randomsoup()
                    else:
                        rulemod = rulemod ^ 2
                        print "rulemod changed to ",rulemod
                        genelife.set_rulemod(rulemod)
                elif event.key == pg.K_s:
                    pg.image.save(screen, "images/genelife_%03d_%08x_%03d.jpeg" % (framenr,repscheme,savecnt))
                    print ("image saved "+"images/genelife_%03d_%08x_%03d.jpeg" % (framenr,repscheme,savecnt))
                    savecnt = savecnt + 1
                elif event.key == pg.K_v:
                    vscrolling=1-vscrolling
                    print "vscrolling changed to ",vscrolling
                    genelife.set_vscrolling()
        if (not mouseclicked and not pause):
            if updatesenabled:
                update_sim(count)
            else:
                colorgrid()
        nspecies=genelife.get_nspecies()
        caption = "Gene Life at step %d coloring %d nspecies %d " % (framenr,colorfunction,nspecies)
        if selection < 8:
            caption = caption + "pairwise selection " + selectiontext0007[selection] + " "
        elif selection<16:
            caption = caption + "LUT encoding " + selectiontext0815[selection-8] + " "
        elif selection<23:
            caption = caption + "multiplane coupling " + selectiontext1623[selection-16] + " "
        if quadrants >= 0:
            paramdat = "repscheme %06x surv. %01x overw. %01x ncoding %06x" % (repscheme,survivalmask,overwritemask,ncoding)
            caption = caption + ("q%1d " % quadrants) + paramdat
        if colorfunction == 4: caption = caption + ("ymax %d " % ymax)
        elif colorfunction == 8: caption = caption + ("offsets (%d,%d,%d) " % (offdx,offdy,offdt))
        if pixeldat: caption = caption + pixeldat
        pg.display.set_caption(caption)
        # pg.transform.scale2x(scr,screen)   # use this for pygame scale2x with smoother
        if scalex2:
            pg.transform.scale2xact(scr,screen)  # use this for custom pygame no smoother
        pg.display.update()                    # copies the screen to the display (or use .flip())
#-----------------------------------------------------------------------------------------------------------

def update_sim(count=True):
    global gol, cgrid
    global golg
    global log2N
    global runparams
    global cnt,framenr,nrun,nskip,ndisp,nhist,nstat

    cnt = cnt+1
    if cnt % ndisp == 0:  # insert the non-displayed iterations & count species
        genelife.genelife_update(nskip, nhist, nstat)
        framenr = framenr + nskip
        if(count):
            genelife.countspecieshash()
    genelife.genelife_update(nrun, nhist, nstat)
    framenr = framenr+nrun
    colorgrid()  # sets  cgrid
    return
#-----------------------------------------------------------------------------------------------------------

def parhelp():
    """ definition of parameters"""
    if selection < 8:
        print "Control bits (left to right):"
        print "____________________________"
        print "Green    ","0. selective birth for 3-live-nbs  ","1. selective birth for 2-live-nbs"
        print "Mid Blue ","2. canonical 0 position vs difft   ","3. bypass selection for 2-live-nbs"
        print "Teal blue","4. enforce birth for 3-live-nbs    ","5. enforce birth for 2-live-nbs"
        print "Red      ","6. 2nd neighbour genetic modulation","7. 1st neighbour genetic masking"
        print "Yellow   ","8. enforce GoL rule if non GoL rule","9. enforce GoL rule last change by non GoL rule */"
        print "Green    ","10-13. allow 2-nb birth only for active subset of 4 canonical configs"
        print "Blue     ","14. Survival for 3-live-nbs        ","15. Survival for 2-live-nbs"
        print "Green    ","16. Gene overwrite for 3-live-nbs  ","17. Gene overwrite for 2-live-nbs"
    elif selection < 12:
        print ""
        if selection < 10:
            print "Control bits for masks to enable gene encoded LUTs (left to right):"
            print "___________________________________________________________________"
            print "Masks for LUT rules for sums 2-6 for survival /left half) and birth (right half)"
            print "blue      ","survival for sum=1-8 separate buttons for the 8 non-zero sums"
            print "green     ","birth    for sum=1-8 separate buttons for the 8 non-zero sums"
        elif selection < 12:
            print "Control bits for masks to enable gene encoded LUTs (left to right):"
            print "___________________________________________________________________"
            print "Distance dept (2 classes) LUT rules for s,se s 2-6 for survival /left half) and birth (right half)"
            print "red      ","sum=2,6 separate buttons for the  3 cases NSEW sum se = 0,1,2 or 2,3,4"
            print "green    ","sum=3,5 separate buttons for the  4 cases NSEW sum se = 0,1,2,3 or 1,2,3,4"
            print "blue     ","sum=4   separate buttons for the 5 cases NSEW sum se = 0,1,2,3,4"
        elif selection < 14:
            print "Control bits for masks to enable gene encoded LUTs (left to right):"
            print "___________________________________________________________________"
            print "A finer division of LUT rules for sums 2-6 for survival /left half) and birth (right half)"
            print "red      ","sum=2,6 separate buttons for the  4 canonical rotations"
            print "green    ","sum=3,5 separate buttons for the  7 canonical rotations"
            print "blue     ","sum=4   separate buttons for the 10 canonical rotations"
        print ""
        print "Additional control bits for selection are in repscheme"
        print "______________________________________________________"
        print "repscheme mod 8 i.e. 0-7 determines selection scheme based on gene"
        print "   # 0 minimum gene as value  # 1 maximum gene as value "
        print "   # 2 minimum number of ones # 3 maximum number of ones"
        print "   # 4 neutral selection # 5 neutral but different selection"
        print "   # 6 penalty function -1 for a survival rule -2 for a birth rule  # 7 not allowed"
        print "repscheme bit 3 (val 0x8) overrides with random choice of ancestor amongst live neighbours"
    print ""
    print "Other controls:"
    print "_______________"
    print "middle mouse","stop simulation [data is retained for possible run() for run/analysis with updatesenabled=True/False]"
    print "left mouse  ","extract information about local state inside the array, or control buttons below"
    print "right mouse ","choose single plane for GoL display in colorfunction 2 for selection 16-19"
    print "<- , ->     ","decrement or increment the colorfunction analysis type mod 9"
    print "h           ","print this help"
    print "<space>     ","pasue simulation"
    print "q,Q         ","incr or decr quadrant parameter choice : -1 = no quadrants, 0-4 are first 5 bit pairs of repscheme, 5,6 surv and overwrite"
    print "r           ","toggle horizon mode on or off: upper half of array obeys unmodified GoL rule"
    print "R           ","toggle random soup domain on or off"
    print "s           ","save current image to file in image subdriectory"
    print "x,X y,Y t,T ","lower (lc) or raise (uc) the (dx,dy,dt) offsets for glider tracking (colorfn 8) (0,0,0)=(all 8 nnb dt=-1)"
    print "v           ","toggle vertical scroll tracking mode : following top most objects and losing lowest objects in contact with 0 row"
    print "+,-         ","increase or decrease ymax for activity display scaled as act/(ymax+act) by a factor of 2"
#-----------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    """ main program with example parameters"""
    nrun=1; ndisp=1000; nskip=0; niter=1;    # simulation time stepping parameters
    nhist = 0                                # set to n to turn on histogram configurations every nth step
    nstat = 0                                # set to n to turn on statistics trace every nth step
    rulemod = runparams[0] = 1               # 0,1 whether to allow GoL rule modifications
                                             # with rulemod 1 2-live-nb birth, 3-live-nb non-birth & non-survival possible
    repscheme = runparams[1] = 8             # repscheme bit 3 (val 0x8) determines whether random choice of ancestor amongst live neighbours
                                             # repscheme mod 8 i.e. 0-7 determines selection scheme based on gene
                                             # 0 minimum gene as value  # 1 maximum gene as value
                                             # 2 minimum number of ones # 3 maximum number of ones
                                             # 4 neutral selection # 5 neutral but different selection
                                             # 6 penalty function -1 for a survival rule -2 for a birth rule  # 7 not allowed
    selection = runparams[2] = 10            # fitness for 2 live neighbor rule : 0-6 see subgenelife.c code
    overwritemask = runparams[3]= 0x00000078 # for selection=10-11 this is the GoL birth mask
    survivalmask = runparams[4] = 0x0000007f # for selection=10-11 this is the GoL survival mask
    colorfunction = runparams[5] = 0         # color function 0(hash), >=1(fnal), 2 nongulstate or color gol planes, 3 notgolrul yellow
                                             # 4 activities 5 genealogy steps 6 genealogy temporal 7 activity scaled colors
    initfield = runparams[6] = 100            # 1 init via 32x32 genepat.dat, n>1 init via nxn rand array
    nlog2pmut = simparams[0] = 8             # log2 gene mutation probability (0 or >56 means no mutation)
    initial1density = simparams[1] =  16384  # initial 1 density in GOL state
                                             # 16384 = nearest to half of guaranteed C rand max value 32767 = 2**15 - 1
    initialrdensity = simparams[2] = 0       # initial density of random genes
    ncoding = simparams[3] = 0               # for selection 10, non zero value means grow plane community from 0
                                             # otherwise (selection<10) no of bits used to encode valid connection functions 1-16
                                             # for selection==8, lut, ncoding 1,2,3 bits per lut entry : 0 implies 3.
    startgenechoice = simparams[4] = 8       # initialize genes to startgene number 0-8 : 8 is random choice of 0-7
    
    genelife.initialize_planes(npoffsets)
    genelife.initialize(runparams,simparams)
    framenr = 0
    cnt=0
    show0()
    # step()
    run()
    
