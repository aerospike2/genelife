// 
// NP 15 Sept 2018:
// including subgenelife.c for low level update routines.
// note 6 params available on command line
// run "actgenelife -h" to see a list of them

// N Packard 18.07.17:
// actgengol.c
// complile:
// cc -o actgengol actgengol.c
// modify output to stdout to output population counts
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

#include <unistd.h>		// for command line parsing


int cmp1 (const void * pa, const void * pb)
{
   // return ( *(int*)pa - *(int*)pb );
   return ((*(uint64_t*)pa > *(uint64_t*)pb)  ? 1 : -1);
}

int cmp2 ( const void *pa, const void *pb )
{
    const uint64_t *a = pa;
    const uint64_t *b = pb;
    if(a[1] == b[1])
        return a[0] > b[0] ? 1 : -1;
    else
        return (int) (b[1] - a[1]);
}


void printspecies(uint64_t golg[]) {  /* counts numbers of all different species using qsort first */
    int ij, k, ijlast, nspecies, counts[N2];
    uint64_t last, golgs[N2];
    uint64_t golgsc[N2][2];
    

    for (ij=0; ij<N2; ij++) { golgs[ij] = golg[ij];  counts[ij] = 0;}  // initialize sorted gene & count arrays to zero

    qsort(golgs, N2, sizeof(uint64_t), cmp1);              // sort in increasing gene order
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
    qsort(golgsc, nspecies, sizeof(golgsc[0]), cmp2);                   // sort in decreasing count order
    for (k=0; k<nspecies; k++) {
        printf("%llx  %llu ",golgsc[k][0],golgsc[k][1]);
    }
    printf("\n");
}

int main (int argc, char *argv[]) {
    int	 i;
    uint64_t *golg;
    int nsteps = 20000;                 // total number of steps to simulate GoL
    
    int myoffs[27] = {0,0,0,
		  -1, 0, 0,
		  -1, 1, 0,
		  0, 1, 0,
		  1, 1, 0,
		  1, 0, 0,
		  1, -1, 0,
		  0, -1, 0,
		  -1, -1, 0};
    int runparams[7];
    int simparams[5];
    int nrunparams=7; int nsimparams=5;

    int opt;
    while ((opt = getopt(argc, argv, "h")) != -1) {
        switch (opt) {
	case 'h':
        default:
	    fprintf(stderr,"Usage: %s  [rulemod=0-1 [repscheme=0-4 [selection=0-2 [nlog2pmut [initial1density [initialrdensity]]]]]\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    state[0] = rand();state[1] = rand();
    Noff = 9;
    runparams[0] = 1;        // 0,1 rulemod
    runparams[1] = 0x19;        // 3-4 repscheme 
    runparams[2] = 1;        // 0-4 selection
    runparams[3] = 3;        // overwritemask = 2 bit mask
    runparams[4] = 0;        // survival mask for two (bit1) and three (bit0) params
    runparams[5] = 1;        // colorfunction 0 or 1
    runparams[6] = 200;        // fileinit 0 or 1, if >1 then size of random init array
    simparams[0] = 8;        // nlog2pmut: gene mutation probability
    simparams[1] = 16384;    // initial1density: nearest to half of guaranteed C rand max value 32767 = 2**15 - 1
    simparams[2] = 0 ;   // initialrdensity: nearest to half of guaranteed C rand max value 32767 = 2**15 - 1
    simparams[3]= 24;        // ncoding
    simparams[4] = 8;        // initialize genes to 0-7 or 8=random

    if (argc>1) runparams[0] = atoi(argv[1]); // if present update rulemod from command line
    if (argc>2) runparams[1] = atoi(argv[2]); // if present update repscheme from command line
    if (argc>3) runparams[2] = atoi(argv[3]); // if present update selection from command line
    if (argc>4) runparams[3] = atoi(argv[4]); // if present update overwritemask from command line
    if (argc>5) runparams[4] = atoi(argv[5]); // if present update survivalmask from command line
    if (argc>6) simparams[0] = atoi(argv[6]); // if present update nlog2pmut from command line
    if (argc>7) simparams[1] = atoi(argv[7]); // if present update initial density from command line
    if (argc>8) simparams[2] = atoi(argv[8]); // if present update init rand density from command line
    if (argc>9) simparams[3] = atoi(argv[9]); // if present update ncoding no of bits used to encode valid connection functions 1-16
    if (argc>10) simparams[4] = atoi(argv[10]); // if present update ncoding no of bits used to encode valid connection functions 1-16
    if (argc>11) runparams[5] = atoi(argv[11]); // if present update colorfunction from command line
    if (argc>12) runparams[6] = atoi(argv[12]); // if present update fileinit from command line

    fprintf(stderr,"Parameters:\n");
    fprintf(stderr,"rulemod-0-1\trepscheme=3,4\tselection=0-4\toverwrite\tsurvival\n");
    fprintf(stderr,"%d\t\t%d\t\t%d\t\t%d\t\t%d\n",runparams[0],runparams[1],runparams[2],runparams[3],runparams[4]);
    fprintf(stderr,"nlog2pmut\tinitial1density\tinitialrdensity\tncoding\tinitgene=0-8\n");
    fprintf(stderr,"%d\t\t%d\t\t%d\t\t%d\t\t%d\n",simparams[0],simparams[1],simparams[2],simparams[3],simparams[4]);

    initialize_planes(myoffs,Noff);
    initialize(runparams,nrunparams,simparams,nsimparams);
    fprintf(stderr,"finished initialize.\n");
    for (i=0; i<nsteps; i++) {                  /* nsteps */
	    golg = planesg[curPlane];
	    printspecies(golg);
	    genelife_update(1,0,0);
    }    
}

