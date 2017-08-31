//
//  subgenelife.c
//  fastgenegol
//
//  Created by John McCaskill on 14.07.17.
//  Copyright © 2017 European Center for Living Technology. All rights reserved.
//
//  The simulation uses the terminal text output as a colour display for the GoL.
//  In order to fit N=128 GoL display on cinema display, use terminal preferences to change font spacings column 1.3 and line 0.65 at 12 pt
//  In order to fit N=128 GoL display on smaller displays, use terminal preferences to change font spacings column 1.0 and line 0.5 at 10 pt
//  To allow escape codes to move cursor, in terminal profiles select "Allow Vt100 application keypad mode"

// From subgenelife2.c
// Modified by Norman Packard Aug 29, 2017
// create ring buffer of planes to accumulate statistics.


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

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

long unsigned int **planes;
long unsigned int **planesg;
int **offsets;
int Noff;
int curPlane = 0;
int newPlane = 1;
int numPlane = 2;
int xL=0,xR=0,yU=0,yD=0;	// offsets for border of stats
long unsigned int *histo;
int numHisto;
                                    // Wikipedia "Xorshift" rewritten here as inline macro &
                                    // Vigna, Sebastiano. "xorshift*/xorshift+ generators and the PRNG shootout". Retrieved 2014-10-25.
#define RAND128P(val) {                                                       \
    long unsigned int x = state[0]; long unsigned int const y = state[1];                       \
	state[0] = y;	x ^= x << 23;  state[1] = x ^ y ^ (x >> 17) ^ (y >> 26);  \
	val = state[1] + y;}

const long unsigned int m1  = 0x5555555555555555; //binary: 0101...           Constants for Hamming distance macro POPCOUNT24C
const long unsigned int m2  = 0x3333333333333333; //binary: 00110011..
const long unsigned int m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const long unsigned int h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...
#define POPCOUNT64C(x, val) {       /* Wikipedia "Hamming Weight" popcount4c alg */  \
    x -= (x >> 1) & m1;             /* put count of each 2 bits into those 2 bits */ \
    x = (x & m2) + ((x >> 2) & m2); /* put count of each 4 bits into those 4 bits */ \
    x = (x + (x >> 4)) & m4;        /* put count of each 8 bits into those 8 bits */ \
    val = (x * h01) >> 56;}         /* left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... */

const long unsigned int hlo = 0xffffffff;
const long unsigned int h02 = h01&hlo;     //lower word multiplier
#define POPCOUNT2X32(x, val1, val2) {      /* Separate Hamming weights of lower and upper 32 bits */  \
    x -= (x >> 1) & m1;             /* put count of each 2 bits into those 2 bits */ \
    x = (x & m2) + ((x >> 2) & m2); /* put count of each 4 bits into those 4 bits */ \
    x = (x + (x >> 4)) & m4;        /* put count of each 8 bits into those 8 bits */ \
    val1 = ((x&hlo) * h02) >> 24;   /* left 8 bits of x + (x<<8) + (x<<16) + (x<<24) for x = x&hlo */ \
    val2 = ((x>>32) * h02) >> 24;}  /* left 8 bits of y + (y<<8) + (y<<16) + (y<<24) for y = x>>32  */

const long unsigned int m3  = 0x0707070707070707; //for cumcount64c selects counter relevant bits only
#define CUMCOUNT64C(x, val) {       /* Assumes gene specifies 8 8-bit counters each with max value 7 */  \
    val = ((x & m3) * h01) >> 56;}  /* left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... */

void countconfigs(){		// count configs specified by offset array
    // each row of the offset array becomes a bit in an address for the histo array.
    int i,j,k,t,x,y,adr,bit;
    unsigned long int *pl;

    for(i=0; i<N; i++){		// rows
	if(i<xL) continue;
	if(i>xR) continue;
	for(j=0; j<N; j++){	// columns
	    if(j<yD) continue;
	    if(j>yU) continue;
	    adr = 0;
	    for(k=0; k<Noff; k++){
		x = j+offsets[k][0];
		y = i+offsets[k][1];
		t = (curPlane+offsets[k][2]) % numPlane;
		pl = planes[t];
		bit = *(pl + y*N +x);
		if(bit>1 || bit <0){
		    fprintf(stderr,"Ack! bit = %d != 0 or 1\n",bit);
		    exit(1);
		}
		adr = (adr<<1) | bit;
	    }
	    histo[adr]++;
	}
    }
}

void get_histo(long unsigned int outhisto[],int numHistoC){
    int i;
    if(numHistoC != numHisto){
	fprintf(stderr,"Ack! numHisto = %d  != numHistoC = %d\n",numHisto,numHistoC);
	exit(1);
    }
    for(i=0; i<numHisto; i++)
	outhisto[i] = histo[i];
}
		    
void init_histo(){
    int i;
    for(i=0; i<numHisto; i++)
	histo[i] = 0;
}

void genelife_update (long unsigned int outgol[], long unsigned int outgolg[], int log2N, int ndsteps, int params[], int N2c, int nparams, int histoflag) {
    /* update GoL for toroidal field which has side length which is a binary power of 2 */
    /* encode without if structures for optimal vector treatment */
    int N = 0x1 << log2N;
    int Nmask = N - 1;            // bit mask for side length, used instead of modulo operation
    int t, k, nmut, nones, nones1, nones2, nones3, nb[8], ij, i, j , jp1, jm1, ip1, im1, p2, p3;
    long unsigned int s, s2or3, nb1i, randnr, randnr1, randnr2, ng, r1, r2, r3, nlog2p, pmask, genediff, birth, newgene;
    long unsigned int genef1,genef2,genef3;
    static float a2=1., a3=1.;
    int nlog2p0 = params[0];
    int nlog2pmut = params[1];
    int selection = params[2];
    int rule2mod = params[3];
    long unsigned int  pmutmask = (0x1 << nlog2pmut) - 1;
    long unsigned int *gol, *newgol, *golg, *newgolg;


    static long unsigned int ngx = 0;
    static int first = 1;
    static long unsigned int state[2]; // State for xorshift pseudorandom number generation. The state must be seeded so that it is not zero

    if (first) {
        state[0] = rand();state[1] = rand();
        first = 0;
    }
    for (t=0; t<ndsteps; t++) {
	gol = planes[curPlane];
	newgol = planes[newPlane];
	golg = planesg[curPlane];
	newgolg = planesg[newPlane];
	for (ij=0; ij<N2c; ij++) {                                               // loop over all sites of 2D torus with side length N
	    i = ij & Nmask;  j = ij >> log2N;                                   // row & column
	    jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
	    ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                         // toroidal i+1, i-1
	    nb[0]=j*N+ip1; nb[1]=j*N+im1; nb[2]=jp1+i; nb[3]=jp1+ip1; nb[4]=jp1+im1; nb[5]=jm1+i; nb[6]=jm1+ip1; nb[7]=jm1+im1; //nbs
	    for (k=0,nb1i=0;k<8;k++) nb1i = (nb1i << (gol[nb[k]]<<2)) + (gol[nb[k]]*k);   // packs non-zero nb indices in first up to 8*4 bits
	    s = gol[nb[0]]+gol[nb[1]]+gol[nb[2]]+gol[nb[3]]+gol[nb[4]]+gol[nb[5]]+gol[nb[6]]+gol[nb[7]]; // number of live nbs
	    s2or3 = (1 - (((s>>3)&1) | ((s>>2)&1))) * (s>>1 & 1);               // 1 if 2 or 3 neighbors are alive
	    if (s2or3 == 1) {                                                   // if 2 or 3 neighbors alive
		RAND128P(randnr);                                                 // expansion inline so compiler recognizes auto-vectorization options
		// compute random neighbor selection ng from one 64-bit random number
		ng = (randnr >> 54) & 0x3;                                          // 0, 1, 2 or 3 with probs each 1/4 : next 4 lines converts this 0,1,2 with prob 1/3
		r3 = (ng == 3 ? 1: 0);                                              // 1 if ng == 3 (invalid value) otherwise zero : prob. is 1/4
		ngx += r3;                                                          // increment external ng counter on such exceptions
		ngx = (ngx == 3 ? 0 : ngx);                                         // modulo 3 counter without division for ng==3 exceptions
		ng = (s&1)?(r3?ngx:ng):(ng&1);                                      // for s==3 use counter value mod 3 if exception, else rand 0,1,2; for s=2 rand 0,1
		newgene = golg[nb[(nb1i>>(ng<<2))& 7]];                             // pick new gene as one of three neighbor
		// compute nones for different selection models
		if (selection == 0) {                                               // neutral model : GoL rule departures depend only on seq diversity
		    genediff = golg[nb[nb1i&0x7]]^golg[nb[(nb1i>>4)&0x7]]^golg[nb[(nb1i>>8)&0x7]];  // gene difference seq based on xor
		    POPCOUNT64C(genediff, nones);}                                    // number of 1s in genediff is Hamming distance
		else if (selection == 1) {                                          // non-neutral model with selection for rule departure probability
		    genediff = newgene;                                               // place the new gene in genediff for count of number of ones
		    POPCOUNT64C(genediff, nones);                                     // number of ones in new gene determines fitness
		    nones = (nones < 16) ? 0 : (nones - 23);}                           // 0 if < 16 otherwise nones-23
		else if (selection == 2) {                                          // non-neutral model based on presence of replicase gene
		    genef1 = golg[nb[nb1i&0x7]];                                      // gene of first neighbor
		    CUMCOUNT64C(genef1, nones);                                       // count of gene determines fitness
		    genef2 = golg[nb[(nb1i>>4)&0x7]];                                 // gene of second neighbor
		    CUMCOUNT64C(genef2, nones1);                                      // number of ones in new gene determines fitness
		    genef3 = golg[nb[(nb1i>>8)&0x7]];                                 // gene of third neighbor
		    CUMCOUNT64C(genef3, nones2);                                      // number of ones in new gene determines fitness
		    nones = nones < nones1 ? nones : nones1;
		    nones = nones < nones2 ? nones : nones2;
		    nones = nones < 16 ? 0 : (nones - 23);}                           // 0 if < 16 otherwise nones-23
		else  {                                                             // full model based on two 32-bit genes for 3nb and 2nb birth
		    genef1 = golg[nb[nb1i&0x7]];                                      // gene of first neighbor
		    POPCOUNT2X32(genef1, nones2, nones3);                             // count of gene ones for 2 32bit words for 3nb and 2nb birth
		    genef2 = golg[nb[(nb1i>>4)&0x7]];                                 // gene of second neighbor
		    POPCOUNT2X32(genef2, nones, nones1);                              // count of gene determines fitness
		    nones2 += nones; nones3 += nones1;
		    genef3 = golg[nb[(nb1i>>8)&0x7]];                                 // gene of third neighbor
		    POPCOUNT2X32(genef3, nones,nones1);                                      // number of ones in new gene determines fitness
		    nones2 += nones; nones3 += nones1;
		    p2=exp(-a2*nones2);
		    p3=exp(-a3*nones3);
                }
		// compute departure and mutation events, mutation position nmut
		nlog2p = nlog2p0 + nones;                                           // real factor alpha not possible here, could do integer conversion
		pmask = (0x1<<nlog2p) - 1L;                                         // probability mask for deviation from gol rules given local hamming
		randnr1 = randnr & pmask;                                           // extract bits from randnr for random trial for 0 on pmask
		r1 = randnr1?0L:1L;                                                   // 1 if lowest nlog2p bits of randnr zero, else zero : i.e. 1 with chance 1/2^nlog2p
		randnr2 = (randnr >> 24) & pmutmask;                                // extract bits from randnr for random trial for 0 on pmask
		r2 = randnr2?0:1;                                                   // 1 if lowest nlog2pmut bits of randnr zero, else zero
		nmut = (randnr >> 48) & 0x3f;                                       // choose mutation position for length 64 gene : from bits 32:37 of randnr
		// complete calculation of newgol and newgolg, including mutation
		newgene = newgene ^ (r2*(0x1L<<nmut));                              // introduce single mutation with probability pmut = probmut
		birth = (0x1L-gol[ij])&((s&1L)^(r1&rule2mod))&0x1;                                 // assuming 2or3 live nbs, birth (value 1) if empty and (s==3 xor r1)
		newgol[ij]  =  gol[ij] | birth ;                                    // new game of life cell value
		newgolg[ij] =  gol[ij]*golg[ij]+birth*newgene;}                     // dies if not 2or3, else old if alive, else new gene if 3 nbs
	    else {                                                              // else not 2 or 3 live neighbors, 0 values
		newgol[ij]  = 0;                                                    // new game of life cell value
		newgolg[ij] = 0;}                                                   // dies if not 2or3
	    emptysites = emptysites + newgol[ij];                               // keep track of empty sites, same information as total activity of occupied sites
	}
	curPlane = (curPlane +1) % numPlane;
	newPlane = (newPlane +1) % numPlane;
	if(histoflag)
	    countconfigs();
	// printf("t = %d nlog2p0=%d nlog2pmut=%d selection=%d rule2mod=%d\n",t,nlog2p0,nlog2pmut,selection,rule2mod);
    }
    for (ij=0; ij<N2c; ij++) {
	outgol[ij] = newgol[ij];        // copy new gol config to output
	outgolg[ij] = newgolg[ij];      // copy new genes to output
    }
}


void printscreen (long unsigned int gol[], long unsigned int golg[], int N, int N2) {   /* print the game of life configuration */
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

void print_gol (long unsigned int gol[], int N, int N2) {   /* print the game of life configuration */
	int	ij;
	for (ij=0; ij<N2; ij++) {
		printf ("%c", gol[ij] ? '*' : ' ');
		if ((ij % N) == N -1) printf ("\n");
	}
}



void initialize_planes(int offs[],  int N) {
    int i,j,idx,Npt;

    if(N%3 !=0)
	fprintf(stderr,"Ouch!  Size of offsets array not a multiple of 3 as expected.");
    Noff = N/3;		// Noff global
    if(Noff>24){
	fprintf(stderr,"Too many offsets!  Max=24");
	exit(1);
    }
    numHisto = 1<<Noff;
    histo = (int *) calloc(numHisto,sizeof(int));

    // install offsets:
    offsets = (int **) calloc(Noff,sizeof(int *));
    for(i=0; i<Noff; i++)
	offsets[i] = (int *) calloc(3,sizeof(int)); // each containing xoff, yoff, toff.
    for(idx=0,i=0; i<Noff; i++){
	for(j=0; j<3; j++){
	    offsets[i][j] = offs[idx];
	    idx++;
	}
    }
    // compute number of planes from toff = 3rd element of each offest vec:
    int tmx = 0;
    int tmn = N;
    int toff,tall;
    for(i=0; i<Noff; i++){
	toff = offsets[i][2];
	if(toff>tmx) tmx = toff;
	if(toff<tmn) tmn = toff;
    }
//    fprintf(stderr,"----------------- txm = %d, tmn = %d",tmx,tmn);
    tall = tmx-tmn;
    numPlane = tall + 2;	// numPlane >= 2

    // compute xR, xL, yU, yD border offsets from offsets matrix
    int mx;
    int mn;
    int off;
    mn=N; mx=0;
    for(i=0; i<Noff; i++){
	off = offsets[i][0];	// X
	if(off>mx) mx = off;
	if(off<mn) mn = off;
    }
    if(mn<0) xL = -mn; else xL=0;
    if(mx>0) xR = N-mx; else xR=N;
    mn=N; mx=0;
    for(i=0; i<Noff; i++){
	off = offsets[i][1];	// Y
	if(off>mx) mx = off;
	if(off<mn) mn = off;
    }
    if(mn<0) yD = -mn; else yD = 0;
    if(mx>0) yU = N-mx; else yU = N;
    

    // initialize planes:
    planes = (long unsigned int **) calloc(numPlane,sizeof(long unsigned int *));
    for(i=0; i<numPlane; i++)
	planes[i] = (long unsigned int *) calloc(N2,sizeof(long unsigned int));
    planesg = (long unsigned int **) calloc(numPlane,sizeof(long unsigned int *));
    for(i=0; i<numPlane; i++)
    planesg[i] = (long unsigned int *) calloc(N2,sizeof(long unsigned int));
    
}

void initialize (int params[], int nparams) {
	int ij;
    int initial1density;
    long unsigned int *gol;
    gol = planes[curPlane];
    static unsigned int rmask = (1 << 15) - 1;

    if (nparams > 4) initial1density = params[4];
    else initial1density = (1 << 14);

    printf("id = %d rmask = %d nparams = %d\n",initial1density, rmask, nparams);
	for (ij=0; ij<N2; ij++) {
		gol[ij] = ((rand() & rmask) < initial1density)?1:0;
	}
}

void initialize_genes (int params[], int nparams) {
    int ij,k;
    long unsigned int g;
    long unsigned int *gol;
    long unsigned int *golg;
    gol = planes[curPlane];
    golg = planesg[curPlane];
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

void countspecies(long unsigned int golg[], int params[], int N2, int nparams) {  /* counts numbers of all different species using qsort first */
    int ij, k, ijlast, nspecies, counts[N2], nones, fitness;
    long unsigned int last, golgs[N2];
    long unsigned int golgsc[N2][2];
    int selection = params[2];
    int nlog2p0 = params[0];

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

