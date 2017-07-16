//
//  fastgol.c
//  fastggenegol
//
//  Created by John McCaskill on 14.07.17.
//  Copyright © 2017 European Center for Living Technology. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

const int log2N = 7;                // toroidal array of side length N = 2 to the power of log2N
const int N = 0x1 << log2N;
const int N2 = N*N;                 // number of sites in toroidal array
const unsigned int Nmask = N - 1;   // bit mask for side length, used instead of modulo operation

const int nlog2p0 = 3;              // p0 = 2 to the power of - nlog2p0
unsigned int p0mask = (0x1 << nlog2p0) - 1;
const int nlog2pmut = 5;            // pmut = probmut = 2 to the power of - nlog2pmut
unsigned int pmutmask = (0x1 << nlog2pmut) - 1;

int nsteps = 10000;
const int ndisp = 10000;

/* The state must be seeded so that it is not zero */
uint64_t state[2];

uint64_t xorshift128plus(void) {
	uint64_t x = state[0];
	uint64_t const y = state[1];
	state[0] = y;
	x ^= x << 23; // a
	state[1] = x ^ y ^ (x >> 17) ^ (y >> 26); // b, c
	return state[1] + y;
}

void update (unsigned int gol[], unsigned int golg[]) {
	/* update GoL for toroidal field which has side length which is a binary power of 2 */
	/* encode without if structures for optimal vector treatment */

	int k, nmut, hamming = 0;
    uint64_t randnr;
    unsigned int s2or3, ij, i, j, jp1, jm1, ip1, im1, s, nb1i, ng, randnr1, randnr2, r1, r2, nlog2p, pmask;
	static unsigned int newgol[N2], newgolg[N2], gene, nb[8];

	for (ij=0; ij<N2; ij++) {        // loop over all sites of 2D torus with side length N
		i = ij & Nmask;              // row
		j = ij >> log2N;             // column
		jp1 = ((j+1) & Nmask)*N;     // toroidal (j+1) * N
		jm1 = ((j-1) & Nmask)*N;     // toroidal (j-1) * N
		ip1 =  (i+1) & Nmask;        // toroidal i+1
		im1 =  (i-1) & Nmask;        // toroidal i-1
        nb[0]=j*N+ip1; nb[1]=j*N+im1; nb[2]=jp1+i; nb[3]=jp1+ip1; nb[4]=jp1+im1; nb[5]=jm1+i; nb[6]=jm1+ip1; nb[7]=jm1+im1; //nbs
        for (k=0,nb1i=0;k<8;k++) nb1i = (nb1i << (gol[nb[k]]<<2)) + (gol[nb[k]]*k); // packs non-zero nb indices in first up to 8*4 bits
		s = gol[nb[0]]+gol[nb[1]]+gol[nb[2]]+gol[nb[3]]+gol[nb[4]]+gol[nb[5]]+gol[nb[6]]+gol[nb[7]]; // number of live nbs
        
        // randnr = xorshift128plus();   //   here we expand inline to allow compiler to optimize auto-vectorization more readily
        uint64_t x = state[0];
        uint64_t const y = state[1];
        state[0] = y;
        x ^= x << 23; // a
        state[1] = x ^ y ^ (x >> 17) ^ (y >> 26); // b, c
        randnr = state[1] + y;
        
        nlog2p = nlog2p0 + hamming;                                         // need to calculate hamming distance
        pmask = p0mask | ((0x1<<hamming) - 1) << nlog2p0;
        randnr1 = randnr & pmask;
        r1 = ((randnr1 | ((randnr1^pmask)+1)) >> nlog2p)&0x1;             // 1 if lowest probpower bits of randnr zero, else zero
        randnr2 = (randnr >> 16) & pmutmask;
        r2 = ((randnr2 | ((randnr2^pmutmask)+1)) >> nlog2pmut)&0x1;         // 1 if lowest probpower bits of randnr zero, else zero
        nmut = (randnr >> 32) & 0x3f;
        ng = ((randnr >> 38) & 0x1)+(s&1)*((randnr >> 39) & 0x1);           // 0, 1 or 2 with probs 1/4,1/2,1/4   : should be 1/3,1/3,1/3
        gene = golg[nb[(nb1i>>(ng<<2))& 7]];                                // pick new gene as one of three neighbor live site genes
        gene = gene ^ (r2*(1<<nmut));
        s2or3 = (1 - (((s>>3)&1) | ((s>>2)&1))) * (s>>1 & 1);               // 1 if 2 or 3 neighbors are alive
        newgol[ij]  = s2or3 * ( (s&1&(1-r1)) | gol[ij] | (r1&(1-(s&1))) );  // new game of life cell value
        newgolg[ij] = s2or3 * (gol[ij]*golg[ij]+(1-gol[ij])*((s&1&(1-r1))|((1-(s&1)&r1)))*gene); // dies if not 2or3, old if alive, else new gene if 3 nbs
	}

	for (ij=0; ij<N2; ij++) {
		gol[ij] = newgol[ij];        // copy new gold config to old grid
        golg[ij] = newgolg[ij];      // copy new genes to old genes
    }
}

void print (unsigned int gol[]) {   /* print the game of life configuration */
	int	ij;
	for (ij=0; ij<N2; ij++) {
		printf ("%c", gol[ij] ? 'x' : ' ');
		if ((ij % N) == N -1) printf ("\n");
	}
}

void initialize (unsigned int gol[]) {
	int ij;
	for (ij=0; ij<N2; ij++) {
		gol[ij] = rand() & 0x1;
	}	
}

void initialize_genes (unsigned int golg[], unsigned int gol[]) {
	int ij,k;
    unsigned int g;
	for (ij=0; ij<N2; ij++) {
        g = 0;
        if (gol[ij] != 0)
            for (k=0; k<32; k++) {
                g = (g << 1) | (rand() & 0x1);
            }
        golg[ij] = g;
        if (golg[ij] == 0 && gol[ij] != 0) printf("zero gene at %d",ij);
	}
    // for (ij=0; ij<20; ij++) printf("gene at %d %u\n",ij,golg[ij]);   // test first 20
}

int cmpfunc (const void * pa, const void * pb)
{
   // return ( *(int*)pa - *(int*)pb );
   return ((*(unsigned int*)pa > *(unsigned int*)pb)  ? 1 : -1);
}

int cmpfunc2 ( const void *pa, const void *pb )
{
    const unsigned int *a = pa;
    const unsigned int *b = pb;
    if(a[1] == b[1])
        return a[0] > b[0] ? 1 : -1;
    else
        return (int) b[1] - (int) a[1];
}

void countspecies(unsigned int golg[]) {  /* counts numbers of all different species using qsort first */
    int ij, k, ijlast, nspecies;
    unsigned int last, golgs[N2], counts[N2];
    unsigned int golgsc[N2][2];
    
    for (ij=0; ij<N2; ij++) { golgs[ij] = golg[ij];  counts[ij] = 0;}  // initialize sorted gene & count arrays to zero

    qsort(golgs, N2, sizeof(unsigned int), cmpfunc);                   // sort in increasing gene order
    for (ij=0,k=0,ijlast=0,last=golgs[0]; ij<N2; ij++) {               // count each new species in sorted list
        if (golgs[ij] != last) {
            last = golgs[ij];
            counts[k++] = ij - ijlast;
            ijlast = ij;
        }
    }
    nspecies = k;  // print excluding 0 since this is most likely an empty site not a true gene
    printf("The number of different non-zero species is %d\n",nspecies);
    
    for (k=0,ij=0;k<nspecies;k++) {     // now condense array to give only different genes with counts
        // printf("species %d with gene %x has counts %d\n",k, golgs[ij],counts[k]);
        golgs[k]=golgs[ij];
        ij = ij + counts[k];
    }
    
    for (k=0; k<nspecies; k++) { golgsc[k][0] = golgs[k];  golgsc[k][1] = counts[k];}  // initialize joint gene & count array
    qsort(golgsc, nspecies, sizeof(golgsc[0]), cmpfunc2);                   // sort in decreasing count order
    for (k=0; k<nspecies; k++) {
        printf("count species %d with gene %x has counts %d\n",k, golgsc[k][0],golgsc[k][1]);
    }
  
}

int main (int argc, char *argv[]) {
    int	 i;
    unsigned int gol[N2], golg[N2];
    
    state[0] = rand();state[1] = rand();
    
    if (argc>1) nsteps = atoi(argv[1]);   /* if present update nsteps from command line */
	initialize (gol);            /* random initial pattern */
    initialize_genes (golg,gol);     /* random initial genes */

    printf("initial pattern  .............................................................................................\n");
    print(gol);
    // for (ij=0; ij<200; ij++) printf("gene at %d %x\n",ij,golg[ij]);   // test first 200
    for (i=0; i<nsteps; i++) {   /* nsteps */
		update (gol, golg);
        if (i%ndisp == 0) {
            printf("after %d steps ...........................................................................................\n",i);
            print (gol);
        }
	}
    printf("and after %d steps ......................................................................................\n",nsteps);
	print (gol);
    printf("and after %d steps ......................................................................................\n",nsteps+1);
    update (gol, golg);                /* one additional step */
    print (gol);
    printf("sort and count of remaining genes ........................................................................\n");
    countspecies(golg);
}
