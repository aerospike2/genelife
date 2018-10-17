// NP 15 Sept 2018:
// including subgenelife.c for low level update routines.
// note 6 params available on command line
// run "actgenelife -h" to see a list of them

// From:
//  fastgol.c
//  fastggenegol
//
//  Created by John McCaskill on 14.07.17.
//  Copyright © 2017 European Center for Living Technology. All rights reserved.
//


#include "subgenelife.c"

#include <unistd.h>

int main (int argc, char *argv[]) {
    int	 i,k;
    uint64_t *gol, *golg;

    int myoffs[27] = {0,0,0,
		  -1, 0, 0,
		  -1, 1, 0,
		  0, 1, 0,
		  1, 1, 0,
		  1, 0, 0,
		  1, -1, 0,
		  0, -1, 0,
		  -1, -1, 0};
    int runparams[5];
    int simparams[5];
    int nrunparams=5; int nsimparams=5;
    int nsteps = 10;	      // total number of steps to simulate GoL
    int ndisp  = 200;	      // display GoL ndisp steps
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

    Noff = 9;
    runparams[0] = 1;          // 0,1 rulemod 
    runparams[1] = 3;        // 0-4 repscheme 
    runparams[2] = 1;        // 0-2 selection 
    runparams[3] = 0;        // NY updated
    runparams[4] = 0;        // NY updated
    simparams[0] = 8;        // nlog2pmut: gene mutation probability
    simparams[1] = 16384;    // initial1density: nearest to half of guaranteed C rand max value 32767 = 2**15 - 1
    simparams[2] = 32768;    // initialrdensity: 32768 makes all genomes active 16384 makes half active, etc.
    simparams[3] = 0;        // NY updated
    simparams[4] = 0;        // NY updated
    
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
    if (argc>11) ndisp = atoi(argv[11]); // if present update init rand density from command line
    if (argc>12) nskip = atoi(argv[12]); // if present update init rand density from command line

    fprintf(stderr,"Parameters:\n");
    fprintf(stderr,"rulemod-0-1\trepscheme=0-4\tselection=0-2\tnlog2pmut\tinitial1density\tinitialrdensity\n");
 fprintf(stderr,"%d\t\t%d\t\t%d\t\t%d\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\n",runparams[0],runparams[1],runparams[2],runparams[3],runparams[4],simparams[0],simparams[1],simparams[2],simparams[3],simparams[4]);
    fprintf(stderr,"ndisp\tnskip\n");
    fprintf(stderr,"%d\t\t%d\n",ndisp,nskip);
    initialize_planes(myoffs,Noff);
    initialize(runparams,nrunparams,simparams,nsimparams);
    for (i=0; i<nsteps; i++) {                  /* nsteps */
	    for(k=0; k< ndisp; k++){
	        gol = planes[curPlane];
	        golg = planesg[curPlane];
	        printxy(gol,golg);
	        genelife_update(1,0);
	    }
	    fprintf(stderr,"finished step %d.\n",i);
	    genelife_update(nskip,0);
	    fprintf(stderr,"step %d\n",i);
    }
}
