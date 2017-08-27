//
//  subgenelife3.c
//  fastgenegol
//
//  Created by John McCaskill and Norman Packard on 11.08.17.
//  Copyright © 2017 European Center for Living Technology. All rights reserved.
//
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define ASCII_ESC 27                // escape for printing terminal commands, such as cursor repositioning

const int log2N = 7;                // toroidal array of side length N = 2 to the power of log2N
const int N = 0x1 << log2N;
const int N2 = N*N;                 // number of sites in toroidal array
const int Nmask = N - 1;            // bit mask for side length, used instead of modulo operation

const int nlog2p0 = 5;              // p0 = 2 to the power of - nlog2p0
const int nlog2pmut = 5;            // pmut = probmut = 2 to the power of - nlog2pmut

int nsteps = 10000;                 // total number of steps to simulate GoL
int ndisp  = 10000;                 // display GoL every ndisp steps
int tdisp  = 0;                     // extra time delay in ms betwene displays
int rule2mod = 1;                   // whether to modify two live nb rule as well or only three nb rule
int selection = 1;                  // fitness model: 0 neutral 1 selected gene prob of rule departure 2 presence of replicase gene
int initial1density = 16384;        // initial density of ones in gol : integer value divide by 2^15 for density
static long unsigned int  emptysites = 0;  // cumulative number of empty sites during simulation updates

const long unsigned int genomemask   = 0xffffffffffffff00;   // 56 bit genome in upper bits
const long unsigned int survgenemask = 0x00000000ffffff00;   // 24 bit gene for survival in lower word bits
const long unsigned int repgenemask  = 0xffffffff00000000;   // 32 bit gene for replication in upper word bits
const long unsigned int sgenomemask  = 0xfffffffffffffff8;   // information passed on in survival
const long unsigned int cntmask      = 0xff;                 // 8 bit counter in lower bits

                                                  // fast routines to count number of ones in various length subwords
const long unsigned int m1  = 0x5555555555555555; //binary: 0101...           Constants for Hamming distance macro POPCOUNTxxC
const long unsigned int m2  = 0x3333333333333333; //binary: 00110011..
const long unsigned int m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const long unsigned int h01 = 0x0101010101010100; //the sum of 256 to the power of 0,1,2,3... for 56 one bit count - both genes
const long unsigned int h02 = 0x0000000001010100; //the sum of 256 to the power of 0,1,2,3... for 24 one bit count - survival
const long unsigned int h03 = 0x0101010100000000; //the sum of 256 to the power of 0,1,2,3... for 32 one bit count - replication
inline void popcount8c(long unsigned x, int *cnt1p) {
    x -= (x >> 1) & m1;             /* put count of each 2 bits into those 2 bits */
    x = (x & m2) + ((x >> 2) & m2); /* put count of each 4 bits into those 4 bits */
    *cnt1p = (x + (x >> 4)) & m4;        /* put count of each 8 bits into those 8 bits */
inline void popcountc(long unsigned x, long unsigned h, int *cnt1p) {
    x -= (x >> 1) & m1;             /* put count of each 2 bits into those 2 bits */
    x = (x & m2) + ((x >> 2) & m2); /* put count of each 4 bits into those 4 bits */
    x = (x + (x >> 4)) & m4;        /* put count of each 8 bits into those 8 bits */
    *cnt1p = (x* h) >> 56;          /* left 8 bits of relevant portion of x + (x<<8) + (x<<16) + (x<<24) + ... */
                                                    // LUT to
static int combi2[256], combi3[256];
void setcombin(int combin[], int n) {
    long unsigned int i,k, cnt;
    for (i=0;k=0,i<256,i++) {
        if (popcount8c(i, &cnt) == n) combin[i] = k++;
        else combin[i] = 64;
    }
}

void genelife_update (long unsigned int gol[], long unsigned int golg[], int log2N, int ndsteps, int params[], int N2c, int nparams) {
	/* update GoL for toroidal field which has side length which is a binary power of 2 */
	/* encode with minimal use of if structures for optimal vector treatment, code optimized with profile */
    int N = 0x1 << log2N;
    int Nmask = N - 1;            // bit mask for side length, used instead of modulo operation
	int t, k, nmut, nones1, nones2, nones3, nb[8], ij, i, j , jp1, jm1, ip1, im1, p2, p3;
    bool s2or3;
    long unsigned int s, nb1i, nbv8, randnr, ng0, ng, ng23, birth, newgene;
    long unsigned int genef1,genef2,genef3;
    static float a2=1., a3=1.;
    int nlog2p0 = params[0];
    int nlog2pmut = params[1];
    int selection = params[2];
    int rule2mod = params[3];
    long unsigned int  pmutmask = (0x1 << nlog2pmut) - 1;
    static long unsigned int  newgol[N2],newgolg[N2];

    static int first = 1;
    static long unsigned int state[2]; // State for xorshift pseudorandom number generation. The state must be seeded so that it is not zero
    
    if (first) {
        state[0] = rand();state[1] = rand();
        first = 0;
    }
    for (t=0; t<ndsteps; t++) {
	  for (ij=0; ij<N2c; ij++) {                                            // loop over all sites of 2D torus with side length N
		i = ij & Nmask;  j = ij >> log2N;                                   // row & column
		jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
		ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                         // toroidal i+1, i-1
        nb[0]=j*N+ip1; nb[1]=j*N+im1; nb[2]=jp1+i; nb[3]=jp1+ip1; nb[4]=jp1+im1; nb[5]=jm1+i; nb[6]=jm1+ip1; nb[7]=jm1+im1; //nbs
        for (k=0,nb1i=0;k<8;k++) nb1i = (nb1i << (gol[nb[k]]<<2)) + (gol[nb[k]]*k);   // packs non-zero nb indices in first up to 8*4 bits
		s = gol[nb[0]]+gol[nb[1]]+gol[nb[2]]+gol[nb[3]]+gol[nb[4]]+gol[nb[5]]+gol[nb[6]]+gol[nb[7]]; // number of live nbs
        if ((s == 2) || (s == 3)) {                                         // if 2 or 3 neighbors alive
            ng0 = golg[ij] & 0x7;                                             // 0-5 value from local rotation counter
            ng = ((ng0+1) == 6)) ? 0 : (ng0+1);                               // increment local rotation counter mod 6
            if (s==2)   ng23 = ng & 0x1;                                      // mod 2 counter for two neighbor case
            else        ng23 = ng >> 1;                                       // mod 3 counter for three neighbor case
            if (gol[ij] == 1) {                                                 // question of survival
                newgene = golg[nb[(nb1i>>(ng23<<2))&0x7]];                        // gene of in phase neighbor
                repage = ((golg[ij] >> 3) + 1) & 0x1f;                            // new survival age of current gene
                popcountc(newgene,h02,&nones1);                                   // ones count of survival gene
                if (nones1 < 3) {                                                 // survives?
                    newgol[ij]  =  gol[ij];                                         // new game of life cell value
                    newgolg[ij] =  (golg[ij]&sgenomemask) | (repage<<3) | ng;}      // genome unchanged but phase & age counters updated
                else {
                    newgol[ij]  =  0UL;
                    newgolg[ij] =  ng;}}
            else {                                                              // question of reproduction
                newgene = golg[nb[(nb1i>>(ng23<<2))&0x7]];                        // gene of in phase neighbor
                popcountc(newgene,h03,&nones1);                                   // ones count of replication gene
                repage = (newgene >> 3) & 0x1f;                                   // age of rep gene
                if (s==3)  {                                                      // certain reproduction for 3 nbs 1
                    if (t&0xf == 0) {                                               // mutation deterministic once every 16 cycles
                        for (k=0,nbv8=0;k<8;k++) nbv8 = (nbv8 << 1) + gol[nb[k]];     // extracts nb values in first 8 bits for mutation
                        newgene = newgene ^ (0x1LU << combi3[nbv8]);                  // mutated gene
                    }
                    newgol[ij]  =  1UL;                                             // new game of life cell value
                    newgolg[ij] =  (newgene & genomemask) | ng;}                    // new gene starts with 0 repage and updated phase ng
                else if (s==2 && (nones1 <= repage) ) {                            // 3 nbs or hamming distance to 0 target less than replication age
                    newgene = 0xffffffffffffff00;
                    newgol[ij]  =  1UL;                                             // new game of life cell value
                    newgolg[ij] =  (newgene & genomemask) | ng;}
                else {
                    newgol[ij]  =  0UL;
                    newgolg[ij] =  ng;}}}
                // still need to increment age
        else {
                newgol[ij]  =  0UL;
                newgolg[ij] =  golg[ij]&0x7;}                                 // preserve only phase counter
        emptysites = emptysites + newgol[ij];                               // keep track of empty sites, same information as total activity of occupied sites
      }

	  for (ij=0; ij<N2c; ij++) {
        gol[ij] = newgol[ij];        // copy new gol config to old one
        golg[ij] = newgolg[ij];      // copy new genes to old genes
      }
    }
}

void initialize (long unsigned int gol[], int params[], int N2, int nparams) {
	int ij;
    int initial1density;
    static unsigned int rmask = (1 << 15) - 1;
    
    if (nparams > 4) initial1density = params[4];
    else initial1density = (1 << 14);
    
    printf("id = %d rmask = %d nparams = %d\n",initial1density, rmask, nparams);
	for (ij=0; ij<N2; ij++) {
		gol[ij] = ((rand() & rmask) < initial1density)?1:0;
	}	
}

void initialize_genes (long unsigned int golg[], long unsigned int gol[], int params[], int N2, int nparams) {
    /* 56 bit gene in upper bits, counters in lower 8 bit */
	int ij,k;
    long unsigned int g;
	for (ij=0; ij<N2; ij++) {
        g = 0;
        if (gol[ij] != 0)
            for (k=0; k<56; k++) {
                g = (g << 1) | (rand() & 0x1);
            }
        golg[ij] = g<<8 + (rand() % 6);    // not optimal uniformity in lowest three bits counter part but this is not critical for initial phase
	}
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

void countspecies(long unsigned int golg[], int params[], int N2, int nparams) {  /* counts numbers of all different species using qsort first */
    int ij, k, ijlast, nspecies, counts[N2], nones, fitness;
    long unsigned int last, golgs[N2];
    long unsigned int golgsc[N2][2];
    int selection = params[2];
    int nlog2p0 = params[0];
    
    for (ij=0; ij<N2; ij++) { golgs[ij] = golg[ij]>>8;  counts[ij] = 0;}  // initialize sorted gene & count arrays, former to gene part, latter to zero

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
        last = golgsc[k][0];
        
        if (selection == 0) {                                               // neutral model : GoL rule departures depend only on seq diversity
                POPCOUNT64C(last, nones);
                fitness = nlog2p0;}
        else if (selection == 1) {                                          // non-neutral model with selection for rule departure probability
                POPCOUNT64C(last, nones);                                     // number of ones in new gene determines fitness
                fitness = nlog2p0 + ((nones < 16) ? 0 : (nones - 23));}       // 0 if < 16 otherwise nones-23
        else {                                                              // non-neutral model based on presence of replicase gene
                POPCOUNT64C(last, nones);                                     // number of ones in new gene determines fitness
                fitness = nlog2p0 + ((nones < 16) ? 0 : (nones - 23));}       // 0 if < 16 otherwise nones-23
        printf("count species %d with gene %lx has counts %lu and %d ones, fitness %d\n",k, golgsc[k][0],golgsc[k][1],nones,fitness);
    }
    printf("cumulative activity = %lu\n",(N2 * (long unsigned int) nsteps) - emptysites);
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
                                    // random numbers not used in this version, but for reference here is our fast PRNG
                                    // Wikipedia "Xorshift" Vigna, Sebastiano. "xorshift*/xorshift+ generators and the PRNG shootout". Retrieved 2014-10-25.
static long unsigned int state[2*N];// State for xorshift pseudorandom number generation. The state must be seeded so that it is not zero
inline long unsigned int rand128p(size_t i) {
    long unsigned int x = state[i]; long unsigned int const y = state[(i+1)];
	state[i] = y;	x ^= x << 23;  state[(i+1)] = x ^ y ^ (x >> 17) ^ (y >> 26);
	return state[(i+1)] + y;}
//  rand128p(&randnr);                                       // use with expansion inline so compiler recognizes auto-vectorization options
