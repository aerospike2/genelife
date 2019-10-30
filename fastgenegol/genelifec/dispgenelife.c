// "dispgenelife.c"
//  project "fastgenegol"
//
//  Created by John McCaskill and Norman Packard 2017-2019.
//  Copyright © 2017 European Center for Living Technology. All rights reserved.
//
//  including subgenelife.c for low level update routines
//  primitive screen graphics in this C calling program (N<=7): for more advanced graphics see python interface
//  note that params are available on command line
//  run "dispgenelife -h" to see a list of them

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#define ASCII_ESC 27                // escape for printing terminal commands, such as cursor repositioning

#include "subgenelife.c"            // set N to 7 for direct screen output in C, significantly larger values possible with python graphics
                                    // see printscreen routine for terminal setup necessary to see direct C output with C calling program
                                    // note that this program was written to run with extended graphical analysis tools from a python notebook

#include <unistd.h>

int main (int argc, char *argv[]) {
    int	 i,k;

    int offsets[24] = {0,0,0, 0,0,-1, 0,0,-2, 0,0,-3, 0,0,-4, 0,0,-5, 0,0,-7};  // offsets used here primarily just to specify the use of 8-plane memory
    int noffsets = 24;
    
    int runparams[10];
    int simparams[6];
    int nrunparams=10;
    int nsimparams=6;
    
    int nsteps = 10;	      // total number of steps to simulate GoL
    int ndisp  = 1000;	      // display GoL ndisp steps
    int nskip = 1000;	      // skip this many

    int opt;
    while ((opt = getopt(argc, argv, "h")) != -1) {
        switch (opt) {
	case 'h':
        default:
	    fprintf(stderr,"Usage: %s  [rulemod=0-1 [repscheme=0-4 [selection=0-2 [nlog2pmut [initial1density [initialrdensity [nstep [nskip]]]]]]]\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }
   
    colorupdate1 = 0;        // do not update color function display or print intermediate results
    runparams[0] = 1;        // 0,1 rulemod
    runparams[1] = 0x170;    // repscheme
    runparams[2] = 8;        // 0-15 selection
    runparams[3] = 0x6;      // overwritemask
    runparams[4] = 0x0;      // survivalmask
    runparams[5] = 0;        // colorfunction 0 to 12, only used in python graphics
    runparams[6] = 0;        // fileinit 0 or 1, if >1 then size of random init array
    runparams[7] = 0x6;      // birthmask
    runparams[8] = 0xff;     // ancselectmask
    runparams[9] = 0;        // colorfunction 2nd window, 0 to 12, only used in python graphics
    
    simparams[0] = 8;        // nlog2pmut: gene mutation probability
    simparams[1] = 8192;     // initial1density: nearest to half of guaranteed C rand max value 32767 = 2**15 - 1
    simparams[2] = 0;        // initialrdensity: 32768 makes all genomes active 16384 makes half active, etc.
    simparams[3] = 1;        // ncodingin: number of coding bits for selection symmetry 8 only 1,2,3 (0 = max i.e. 3)
    simparams[4] = 8;        // startgenechoice (8 for random choice between all 8 possible)
    simparams[5] = 1234567;  // ranseed initialization
    
    if (argc>1) runparams[0] = atoi(argv[1]); // if present update rulemod from command line
    if (argc>2) runparams[1] = atoi(argv[2]); // if present update repscheme from command line
    if (argc>3) runparams[2] = atoi(argv[3]); // if present update selection from command line
    if (argc>4) runparams[3] = atoi(argv[4]); // if present update selection from command line
    if (argc>5) runparams[4] = atoi(argv[5]); // if present update selection from command line
    if (argc>6) simparams[0] = atoi(argv[6]); // if present update nlog2pmut from command line
    if (argc>7) simparams[1] = atoi(argv[7]); // if present update initial density from command line
    if (argc>8) simparams[2] = atoi(argv[8]); // if present update init rand density from command line
    if (argc>9) simparams[3] = atoi(argv[9]); // if present update init rand density from command line
    if (argc>10) simparams[4] = atoi(argv[10]); // if present update init rand density from command line
    if (argc>11) runparams[5] = atoi(argv[11]); // if present update colorfunction from command line
    if (argc>12) runparams[6] = atoi(argv[12]); // if present update fileinit from command line
    if (argc>13) ndisp = atoi(argv[13]); // if present update init rand density from command line
    if (argc>14) nskip = atoi(argv[14]); // if present update init rand density from command line

    fprintf(stderr,"Parameters:\n");
    fprintf(stderr,"rulemod-0-1\trepscheme=0-4\tselection=0-2\tnlog2pmut\tinitial1density\tinitialrdensity\n");
    fprintf(stderr,"%d\t\t%d\t\t%d\t\t%d\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\n",runparams[0],runparams[1],runparams[2],runparams[3],runparams[4],
         simparams[0],simparams[1],simparams[2],simparams[3],simparams[4]);
    fprintf(stderr,"ndisp\tnskip\n");
    fprintf(stderr,"%d\t\t%d\n",ndisp,nskip);
    initialize_planes(offsets,noffsets);
    initialize(runparams,nrunparams,simparams,nsimparams);
    for (i=0; i<nsteps; i++) {                  /* nsteps */
	    for(k=0; k< ndisp; k++){
            // printxy(gol,golg); // simple character output
	        printscreen(gol,golg);  // colour xterm output moving cursor with esc codes
	        genelife_update(1,0,0);
            usleep(20000);
	    }
	    fprintf(stderr,"finished step %d.\n",i);
	    genelife_update(nskip,0,0);
	    fprintf(stderr,"step %d\n",i);
    }
}
