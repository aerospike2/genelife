// 
// N Packard 18.07.17:
// actgengol.c
// complile:
// cc -o actgengol actgengol.c
// modify output to stdout to ouput population counts
// to run with activity.py:
// activity.py actgengol
//
// From:
//  fastgol.c
//  fastggenegol
//
//  Created by John McCaskill on 14.07.17.
//  Copyright © 2017 European Center for Living Technology. All rights reserved.
//

#include "subgenelife.c"

void printspecies(long unsigned int golg[]) {  /* counts numbers of all different species using qsort first */
    int ij, k, ijlast, nspecies, counts[N2];
    long unsigned int last, golgs[N2];
    long unsigned int golgsc[N2][2];
    
    for (ij=0; ij<N2; ij++) { golgs[ij] = golg[ij];  counts[ij] = 0;}  // initialize sorted gene & count arrays to zero

    qsort(golgs, N2, sizeof(long unsigned int), cmpfunc);              // sort in increasing gene order
    for (ij=0,k=0,ijlast=0,last=golgs[0]; ij<N2; ij++) {               // count each new species in sorted list
        if (golgs[ij] != last) {
            last = golgs[ij];
            counts[k++] = ij - ijlast;
            ijlast = ij;
        }
    }
    nspecies = k;  // print excluding 0 since this is most likely an empty site not a true gene
    for (k=0,ij=0;k<nspecies;k++) {     // now condense array to give only different genes with counts
        // printf("species %4d with gene %x has counts %d\n",k, golgs[ij],counts[k]);
        golgs[k]=golgs[ij];
        ij = ij + counts[k];
    }
    
    for (k=0; k<nspecies; k++) { golgsc[k][0] = golgs[k];  golgsc[k][1] = counts[k];}  // initialize joint gene & count array
    qsort(golgsc, nspecies, sizeof(golgsc[0]), cmpfunc2);                   // sort in decreasing count order
    for (k=0; k<nspecies; k++) {
        printf("%lx  %lu ",golgsc[k][0],golgsc[k][1]);
    }
    printf("\n");
}

int main (int argc, char *argv[]) {
    int	 i,k;
    long unsigned int *golg;
    int nsteps = 10000;                 // total number of steps to simulate GoL
    
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

    state[0] = rand();state[1] = rand();
    Noff = 9;
    runparams[0] = 1;          // 0,1 rulemod 
    runparams[1] = 3;        // 0-4 repscheme 
    runparams[2] = 1;        // 0-2 selection 

    simparams[0] = 8;        // nlog2p0:   base prob of GOL departure 1/2^nlog2p0 
    simparams[1] = 8;        // nlog2pmut: gene mutation probability 
    simparams[2] = 16384;   // initial1density: nearest to half of guaranteed C rand max value 32767 = 2**15 - 1 

    if (argc>1) simparams[0] = atoi(argv[1]);         /* if present update nlog2p0 from command line */
    if (argc>2) simparams[1] = atoi(argv[2]);          /* if present update nlog2pmut from command line */
    if (argc>3) simparams[2] = atoi(argv[3]);          /* if present update initialdensity from command line */

    initialize_planes(myoffs,Noff);
    initialize(runparams,nrunparams,simparams,nsimparams);
    initialize_genes(simparams,nsimparams);
    fprintf(stderr,"finished initialize.\n");
    for (i=0; i<nsteps; i++) {                  /* nsteps */
	golg = planesg[curPlane];
	printspecies(golg);
	genelife_update(1,0);
    }    
}

