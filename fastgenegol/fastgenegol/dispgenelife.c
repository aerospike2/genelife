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
    long unsigned int *gol;
    long unsigned int *golg;

    int myoffs[27] = {0,0,0,
		  -1, 0, 0,
		  -1, 1, 0,
		  0, 1, 0,
		  1, 1, 0,
		  1, 0, 0,
		  1, -1, 0,
		  0, -1, 0,
		  -1, -1, 0};
    int runparams[3];
    int simparams[3];
    int nrunparams=3; int nsimparams=3;
    int nsteps = 10;	      // total number of steps to simulate GoL
    int ndisp  = 100;	      // display GoL ndisp steps
    int nskip = 1000;	      // skip this many

    int opt;
    while ((opt = getopt(argc, argv, "h")) != -1) {
        switch (opt) {
	case 'h':
        default:
	    fprintf(stderr,"Usage: %s  [rulemod=0-1 [repscheme=0-4 [selection=0-2 [nlog2p0 [nlog2pmut [initialdensity]]]]]]\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    char strh[10];

    // strcpy(strh,"-h");
    // if(argc>0){
    // 	if(strcmp(strh,argv[1])==0){
    // 	    
    // 	    exit(0);
    // 	}
    // }



    Noff = 9;
    runparams[0] = 1;          // 0,1 rulemod 
    runparams[1] = 3;        // 0-4 repscheme 
    runparams[2] = 1;        // 0-2 selection 

    simparams[0] = 8;        // nlog2p0:   base prob of GOL departure 1/2^nlog2p0 
    simparams[1] = 8;        // nlog2pmut: gene mutation probability 
    simparams[2] = 16384;   // initial1density: nearest to half of guaranteed C rand max value 32767 = 2**15 - 1 

    if (argc>1) runparams[0] = atoi(argv[1]); // if present update rulemod from command line 
    if (argc>2) runparams[1] = atoi(argv[2]); // if present update repscheme from command line 
    if (argc>3) runparams[2] = atoi(argv[3]); // if present update selection from command line 
    if (argc>4) simparams[0] = atoi(argv[4]); // if present update nlog2p0 from command line 
    if (argc>5) simparams[1] = atoi(argv[5]); // if present update nlog2pmut from command line 
    if (argc>6) simparams[2] = atoi(argv[6]); // if present update initialdensity from command line 

    initialize_planes(myoffs,Noff);
    initialize(runparams,nrunparams,simparams,nsimparams);
    initialize_genes(simparams,nsimparams);
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
