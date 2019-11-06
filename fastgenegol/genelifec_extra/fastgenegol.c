//
//  fastgenegol.c
//  fastgenegol
//
//  Created by John McCaskill on 14.07.17.
//  Copyright © 2017 European Center for Living Technology. All rights reserved.
//
//  The simulation uses the terminal text output as a colour display for the GoL.
//  In order to fit N=128 GoL display on cinema display, use terminal preferences to change font spacings column 1.3 and line 0.65 at 12 pt
//  In order to fit N=128 GoL display on smaller displays, use terminal preferences to change font spacings column 1.0 and line 0.5 at 10 pt
//  To allow escape codes to move cursor, in terminal profiles select "Allow Vt100 application keypad mode"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#define ASCII_ESC 27                // escape for printing terminal commands, such as cursor repositioning

const int log2N = 7;                // toroidal array of side length N = 2 to the power of log2N
const int N = 0x1 << log2N;
const int N2 = N*N;                 // number of sites in toroidal array
const int Nmask = N - 1;            // bit mask for side length, used instead of modulo operation

const int nlog2p0 = 8;              // p0 = 2 to the power of - nlog2p0
const int nlog2pmut = 5;            // pmut = probmut = 2 to the power of - nlog2pmut

int nsteps = 10000;                 // total number of steps to simulate GoL
int ndisp = 10000;                  // display GoL every ndisp steps
int tdisp = 0;                      // extra time delay in ms betwene displays
static long unsigned int  emptysites = 0;  // cumulative number of empty sites during simulation updates

long unsigned int state[2];                  // State for xorshift pseudorandom number generation. The state must be seeded so that it is not zero
                                    // Wikipedia "Xorshift" rewritten here as inline macro &
                                    // Vigna, Sebastiano. "xorshift*/xorshift+ generators and the PRNG shootout". Retrieved 2014-10-25.
#define RAND128P(val) {                                                       \
    long unsigned int x = state[0]; long unsigned int const y = state[1];                       \
	state[0] = y;	x ^= x << 23;  state[1] = x ^ y ^ (x >> 17) ^ (y >> 26);  \
	val = state[1] + y;}

const long unsigned int m1  = 0x5555555555555555; //binary: 0101...           Constants for Hamming distance macro POPCOUNT64C
const long unsigned int m2  = 0x3333333333333333; //binary: 00110011..
const long unsigned int m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const long unsigned int h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

#define POPCOUNT64C(x, val) {          /* Wikipedia "Hamming Weight" popcount4c alg */  \
    gx = x;                            /* since the macro modifies gx, make a copy first */\
    gx -= (gx >> 1) & m1;              /* put count of each 2 bits into those 2 bits */ \
    gx = (gx & m2) + ((gx >> 2) & m2); /* put count of each 4 bits into those 4 bits */ \
    gx = (gx + (gx >> 4)) & m4;        /* put count of each 8 bits into those 8 bits */ \
    val = (gx * h01) >> 56;}           /* left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... */

void update (long unsigned int gol[], long unsigned int golg[]) {
	/* update GoL for toroidal field which has side length which is a binary power of 2 */
	/* encode without if structures for optimal vector treatment */

	int k, nmut, hamming, nb[8], ij, i, j , jp1, jm1, ip1, im1;
    long unsigned int s, s2or3, nb1i, randnr, randnr1, randnr2, ng, r1, r2, r3, nlog2p, pmask, gx, genediff, birth, newgene;
    static long unsigned int  newgol[N2],newgolg[N2];
    static long unsigned int  pmutmask = (0x1 << nlog2pmut) - 1;
    static long unsigned int ngx = 0;


	for (ij=0; ij<N2; ij++) {                                               // loop over all sites of 2D torus with side length N
		i = ij & Nmask;  j = ij >> log2N;                                   // row & column
		jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
		ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                         // toroidal i+1, i-1
        nb[0]=j*N+ip1; nb[1]=j*N+im1; nb[2]=jp1+i; nb[3]=jp1+ip1; nb[4]=jp1+im1; nb[5]=jm1+i; nb[6]=jm1+ip1; nb[7]=jm1+im1; //nbs
        for (k=0,nb1i=0;k<8;k++) nb1i = (nb1i << (gol[nb[k]]<<2)) + (gol[nb[k]]*k);   // packs non-zero nb indices in first up to 8*4 bits
		s = gol[nb[0]]+gol[nb[1]]+gol[nb[2]]+gol[nb[3]]+gol[nb[4]]+gol[nb[5]]+gol[nb[6]]+gol[nb[7]]; // number of live nbs
        s2or3 = (1 - (((s>>3)&1) | ((s>>2)&1))) * (s>>1 & 1);               // 1 if 2 or 3 neighbors are alive
        genediff = s2or3 * (golg[nb[nb1i&0x7]]^golg[nb[(nb1i>>4)&0x7]]^golg[nb[(nb1i>>8)&0x7]]);  // gene difference sequence based on xor for 2 or 3 live nbs
        POPCOUNT64C(genediff, hamming);
        
        // compute random events r1 and r2, neighbor selection ng and single mutation position nmut from same 64-bit random number
        RAND128P(randnr);                                                   // expansion inline so compiler recognizes auto-vectorization options
        nlog2p = nlog2p0 + hamming;                                         // real factor alpha not possible here, could do integer conversion
        pmask = (0x1<<nlog2p) - 1;                                          // probability mask for deviation from gol rules given local hamming
        randnr1 = randnr & pmask;                                           // extract bits from randnr for random trial for 0 on pmask
        r1 = ((randnr1 | ((randnr1^pmask)+1)) >> nlog2p)&0x1;               // 1 if lowest nlog2p bits of randnr zero, else zero : i.e. 1 with chance 1/2^nlog2p
        randnr2 = (randnr >> 16) & pmutmask;                                // extract bits from randnr for random trial for 0 on pmask
        r2 = ((randnr2 | ((randnr2^pmutmask)+1)) >> nlog2pmut)&0x1;         // 1 if lowest nlog2pmut bits of randnr zero, else zero
        nmut = (randnr >> 32) & 0x3f;                                       // choose mutation position for length 64 gene
        ng = (randnr >> 38) & 0x3;                                          // 0, 1, 2 or 3 with probs each 1/4 : next 4 lines converts this 0,1,2 with prob 1/3
        r3 = ((ng+1) >> 2) & 0x1;                                           // 1 if ng == 3 (invalid value) otherwise zero : prob. is 1/4
        ngx += r3;                                                          // increment external ng counter on such exceptions
        ngx = (1-(((ngx+1)>> 2)&0x1))*ngx;                                  // modulo 3 counter without division for ng==3 exceptions
        ng = (s&1)*((1-r3)*ng+r3*ngx)+(1-(s&1))*(ng&1);                     // use counter value mod 3 if exception, else random 0,1,2 for s==3; random 0,1 for 2
        
                                  // number of 1s in genediff is Hamming distance

        newgene = golg[nb[(nb1i>>(ng<<2))& 7]];                             // pick new gene as one of three neighbor live site genes
        newgene = newgene ^ (r2*(0x1L<<nmut));                              // introduce single mutation with probability pmut = probmut
        birth = (1-gol[ij])&((s&1)^r1)&0x1;                                 // assuming 2or3 live nbs, birth (value 1) if empty and (s==3 xor r1)
        newgol[ij]  = s2or3 * ( gol[ij] | birth );                          // new game of life cell value
        newgolg[ij] = s2or3 * ( gol[ij]*golg[ij]+birth*newgene);            // dies if not 2or3, else old if alive, else new gene if 3 nbs

        emptysites = emptysites + newgol[ij];
    }

	for (ij=0; ij<N2; ij++) {
        gol[ij] = newgol[ij];        // copy new gol config to old one
        golg[ij] = newgolg[ij];      // copy new genes to old genes
    }
}

void printscreen (long unsigned int gol[], long unsigned int golg[]) {   /* print the game of life configuration */
	int	ij, col;
    // https://stackoverflow.com/questions/27159322/rgb-values-of-the-colors-in-the-ansi-extended-colors-index-17-255
    printf("\e[38;5;255;48;5;238m");
	for (ij=0; ij<N2; ij++) {
        col = 32+((golg[ij]>>57)&0x7f);
		printf ("\e[38;5;%dm%c", col, gol[ij] ? '*' : ' ');
		if ((ij % N) == N -1) printf ("\n");
	}
    printf("\e[38;5;238;48;5;255m");
}

void print (long unsigned int gol[], long unsigned int golg[]) {   /* print the game of life configuration */
	int	ij;
	for (ij=0; ij<N2; ij++) {
		printf ("%c", gol[ij] ? '*' : ' ');
		if ((ij % N) == N -1) printf ("\n");
	}
}

void initialize (long unsigned int gol[]) {
	int ij;
	for (ij=0; ij<N2; ij++) {
		gol[ij] = rand() & 0x1;
	}	
}

void initialize_genes (long unsigned int golg[], long unsigned int gol[]) {
	int ij,k;
    uint64_t g;
	for (ij=0; ij<N2; ij++) {
        g = 0;
        if (gol[ij] != 0)
            for (k=0; k<64; k++) {
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
   return ((*(long unsigned int*)pa > *(long unsigned int*)pb)  ? 1 : -1);
}

int cmpfunc2 ( const void *pa, const void *pb )
{
    const long unsigned int *a = pa;
    const long unsigned int *b = pb;
    if(a[1] == b[1])
        return a[0] > b[0] ? 1 : -1;
    else
        return (int) (b[1] - a[1]);
}

void countspecies(long unsigned int golg[]) {  /* counts numbers of all different species using qsort first */
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
    printf("The number of different non-zero species is %d\n",nspecies);
    
    for (k=0,ij=0;k<nspecies;k++) {     // now condense array to give only different genes with counts
        // printf("species %4d with gene %x has counts %d\n",k, golgs[ij],counts[k]);
        golgs[k]=golgs[ij];
        ij = ij + counts[k];
    }
    
    for (k=0; k<nspecies; k++) { golgsc[k][0] = golgs[k];  golgsc[k][1] = counts[k];}  // initialize joint gene & count array
    qsort(golgsc, nspecies, sizeof(golgsc[0]), cmpfunc2);                   // sort in decreasing count order
    for (k=0; k<nspecies; k++) {
        printf("count species %d with gene %llx has counts %llu\n",k, golgsc[k][0],golgsc[k][1]);
    }
    printf("cumulative activity = %llu\n",(N2 * (uint64_t) nsteps) - emptysites);
}

void delay(int milliseconds)
{
    long pause;
    clock_t now,then;
    
    pause = milliseconds*(CLOCKS_PER_SEC/1000);
    now = then = clock();
    while( (now-then) < pause )
        now = clock();
}

int main (int argc, char *argv[]) {
    int	 i;
    long unsigned int gol[N2];
    long unsigned int golg[N2];
    
    state[0] = rand();state[1] = rand();
    
    if (argc>1) nsteps = atoi(argv[1]);         /* if present update nsteps from command line */
    if (argc>2) ndisp = atoi(argv[2]);          /* if present update ndisp from command line */
    if (argc>3) tdisp = atoi(argv[3]);          /* if present update tdisp from command line */

	initialize (gol);                           /* random initial pattern */
    initialize_genes (golg,gol);                /* random initial genes */

    printf("initial pattern  ..............................................................................................................\n");
    print(gol,golg);
    for (i=0; i<nsteps; i++) {                  /* nsteps */
		update (gol, golg);                     /* update game of life and gene array by one step */
        if ((ndisp < nsteps) && (i%ndisp == 0)) {                     /* only display as required */
            printf("%d steps\n",i);             /* display step counter */
            printscreen (gol,golg);                   /* print genetic game of life using colour coded genes for live states */
            printf( "%c[%dA", ASCII_ESC, N+1 ); /* move cursor up N+1 lines on VT100 screen */
            if (tdisp) delay(tdisp);            /* delay to allow display to keep up if desired */
        }
	}
    printf("and after %d steps ............................................................................................................\n",nsteps);
	print (gol,golg);
    printf("sort and count of remaining genes ..............................................................................................\n");
    countspecies(golg);
}
