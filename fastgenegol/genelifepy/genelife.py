import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpcolors
from matplotlib.colors import ListedColormap
import matplotlib.animation as animation
import matplotlib
import genelife_update_module as genelife


log2N = 8                       #  log2N=7 => N=128
N = 2**log2N
N2 = N*N
Nmask = N-1
gol = np.zeros(N2,np.uint64)
golg = np.zeros(N2,np.uint64)
cgrid = np.zeros((N,N),np.uint)



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

my_cmap = rand_cmap(257, type='bright', first_color_black=True, last_color_black=False, verbose=False)

def colorgrid(N):
    """ colors array according to grid and genegrid using colormethod"""
    global gol,cgrid,golg
    for i in range(N):                 # make this simple array copy operation more efficient
        for j in range(N):
            if gol[i+j*N]:
                cgrid[i,j] = 1 + (int(golg[i+j*N]//0x100000000000000))
            else:
                cgrid[i,j] = 0
    return



cnt = 0
framenr = 0

# doanimation configured for:
#     |--ndisp--|--------nskip-------|  repeated niter times


mat = []
ndisp = 100
def doanimation(nrun=1,         # number of CA iterations per animation time step
                ndisp=100,      # number of display steps
                nskip = 500,    # CA iterations to skip between display
                niter=8):       # number of ndisp-nskip steps
    # use my_dpi for your monitor's dpi, to get cells to line up on pixel boundaries.
    # for my_dpi go to http://www.infobyip.com/detectmonitordpi.php
    my_dpi=96
    #fig, ax = plt.subplots(figsize=(10,10))
    fig, ax = plt.subplots(figsize=((3*256)/my_dpi, (3*256)/my_dpi), dpi=my_dpi, frameon=False)
    ax.axis('off')
    cnt=0
    framenr = 0
    time_text = ax.text(0.05, 0.95,'',horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)
    colorgrid(N)
    def update_anim(data):
        global gol, cgrid
        global cnt
        global framenr
        cnt = cnt+1
        if cnt % ndisp == 0:  # insert the non-displayed iterations
            genelife.genelife_update(nskip, 0)
        genelife.genelife_update(1, 0) # 1 should be nrun
        genelife.get_curgol(gol)
        genelife.get_curgolg(golg)
        framenr = framenr+1         # 1 should be nrun
        colorgrid(N)                # sets cgrid from gol, golg
        mat.set_data(cgrid)
        time_text.set_text('cnt = %.1d' % framenr)
        time_text.set_color('w')

    mat = ax.matshow(cgrid, cmap=my_cmap, vmin=0.01, vmax=257)  # was vmax = LEN+1
    ani = animation.FuncAnimation(fig, update_anim, interval=1, # ~40 ms to see every time step
                                      save_count=0, frames=niter*ndisp, repeat = False)
    plt.show()

if __name__ == '__main__':
    offsets =  [[ 0, 0, 0],
                [-1, 0, 0],
                [-1, 1, 0],
                [ 0, 1, 0],
                [ 1, 1, 0],
                [ 1, 0, 0],
                [ 1,-1, 0],
                [ 0,-1, 0],
                [-1,-1, 0]]
    cgrid = np.zeros((N,N),np.uint)

    runparams = np.zeros(3,np.int32)    # 3 parameters passed to C
    rulemod = runparams[0] = 1          # 0,1
    repscheme = runparams[1] = 4        # 0-4
    selection = runparams[2] = 1        # 0-2
    
    simparams = np.zeros(5,np.int32)    # 5 parameters passed to C
    nlog2p0   = simparams[0] = 8        # base prob of GOL departure 1/2^nlog2p0
    nlog2pmut = simparams[1] = 8        # gene mutation probability
    nloglog2p1 = simparams[2] = 2       # power of 2 reduction in influence of nr of 1s in seq on prob
    initial1density = simparams[3] = 16384   # nearest to half of guaranteed C rand max value 32767 = 2**15 - 1
    initialrdensity = simparams[4] = 16384   # nearest to half of guaranteed C rand max value 32767 = 2**15 - 1
    flatoff =  [x for sublist in offsets for x in sublist]
    npoffsets = np.array(flatoff,np.int32)
    genelife.initialize_planes(npoffsets)
    genelife.initialize(runparams,simparams)
    genelife.initialize_genes(runparams)              // correction by John, but mistake had no effect

    doanimation(nrun=1,ndisp=100,nskip=1000,niter=10)

#    genelife.get_curgolg(golg)
#    genelife.countspecies(golg, simparams)
    
