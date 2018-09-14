#include "subgenelife.c"

/**
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>


extern long unsigned int **planes;
extern long unsigned int **planesg;
extern int curPlane;
extern int Noff;
extern void initialize_planes(int *, int);
extern void initialize(int *, int, int*, int);
extern void initialize_genes(int *, int);
extern void printxy(long unsigned int *,long unsigned int *);
extern void genelife_update(int,int);
**/

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
