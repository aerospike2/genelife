// Subsequent fix of new integer types and error in {} structure John McCaskill, Sep 4, 2017
// 
// From subgenelife.c
// Modified by John McCaskill for coded departures from GoL Sep 21, 2018
//
//  subgenelife_codedep.c
//  fastgenegol
//
//  Created by John McCaskill on 14.07.17.
//  Copyright © 2017 European Center for Living Technology. All rights reserved.
//



#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#define ASCII_ESC 27                // escape for printing terminal commands, such as cursor repositioning

const int log2N = 8;                // toroidal array of side length N = 2 to the power of log2N
const int N = 0x1 << log2N;
const int N2 = N*N;                 // number of sites in toroidal array
const int Nmask = N - 1;            // bit mask for side length, used instead of modulo operation

// from fastgenelifed.c (9/9/18):
                                    // integer negative log 2 of max prob p0 for rule departure and prob pmut of mutation
                                    // the maximum prob is further reduced by a sequence-specific prob p1 (p=p0*p1)
                                    // p1 goes from 1 down to 2^-32 depending on genetic differences
int nlog2pmut = 5;                  // pmut (prob of mutation) = 2^(-nlog2pmut), minimum
long unsigned int pmutmask;         // binary mask so that prob of choosing zero is pmut, assigned

//const int nlog2pmut = 5;            // pmut = probmut = 2 to the power of - nlog2pmut

int totsteps=0;
int rulemod = 1;                    // det: whether to modify GoL rule for 2 and 3 live neighbours : opposite outcome with small probability p0
int selection = 1;                  // fitness model: 0 neutral 1 selected gene prob of rule departure 2 presence of replicase gene
int repscheme = 1;                  // replication scheme: 0 random choice , 1 XOR of all 3, 2 consensus, 3 unique det choice, 4 most different

int initial1density = (1<<15)>>1;        // initial density of ones in gol as integer value, divide by 2^15 for true density
int initialrdensity = (1<<15)>>1;        // initial density of random genes in live sites, divide by 2^15 for true density
static long unsigned int  emptysites = 0;  // cumulative number of empty sites during simulation updates

int Noff = 9;                           // number of offsets
int **offsets;                      // array of offsets (2D + time) for planes
int curPlane = 0;
int newPlane = 1;
int xL=0,xR=0,yU=0,yD=0;            // offsets for border of stats
long unsigned int *histo;
int numHisto;

// initialize planes:
int numPlane = 8;
long unsigned int *planes[8];         // ring buffer planes of gol array states
long unsigned int *planesg[8];        // ring buffer planes of golg genes
long unsigned int plane0[N2];
long unsigned int plane1[N2];
long unsigned int plane2[N2];
long unsigned int plane3[N2];
long unsigned int plane4[N2];
long unsigned int plane5[N2];
long unsigned int plane6[N2];
long unsigned int plane7[N2];
long unsigned int planeg0[N2];
long unsigned int planeg1[N2];
long unsigned int planeg2[N2];
long unsigned int planeg3[N2];
long unsigned int planeg4[N2];
long unsigned int planeg5[N2];
long unsigned int planeg6[N2];
long unsigned int planeg7[N2];

                                    // Wikipedia "Xorshift" rewritten here as inline macro &
                                    // Vigna, Sebastiano. "xorshift*/xorshift+ generators and the PRNG shootout". Retrieved 2014-10-25.
static long unsigned int state[2]; // State for xorshift pseudorandom number generation. The state must be seeded so that it is not zero
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
    int i,j,k,t,x,y;
    unsigned long int *pl, adr, bit;

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
		        if(bit!=1 || bit != 0){                                            // J nmodified check for unsigned values */
		            fprintf(stderr,"Ack! bit = %lud != 0 or 1\n",bit);
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
    for(i=0; i<numHisto; i++) outhisto[i] = histo[i];
}

void init_histo(){     // initialize the history array to zero
    int i;
    for(i=0; i<numHisto; i++)        histo[i] = 0;
}

void update(long unsigned int gol[], long unsigned int golg[],long unsigned int newgol[], long unsigned int newgolg[]){
    /* update GoL for toroidal field which has side length which is a binary power of 2 */
    /* encode without if structures for optimal vector treatment */
    /* use only with repscheme == 2 until further notice, later repscheme == 4 and others will be supported */

    int k, kanc, nmut, l;
    int nb[8], ij, i, j, jp1, jm1, ip1, im1;
    int nb1[8], ij1, i1, j1, j1p1, j1m1, i1p1, i1m1;
    long unsigned int s, gs, s2or3, s2, sl, nbi, nb1i, nbil, randnr, randnr2, r2;
    long unsigned int rulemodl;
    long unsigned int gene, newgene, genematch, genelink;
    // long unsigned int andgene, orgene;
    // int nones;
    long unsigned int survive, birth;
    unsigned int genebyte, cmask;
    static long unsigned matching = 1;

    totsteps++;
    randnr = 0x0123456789abcdef;                                            // initialize random nr
    rulemodl = rulemod ? 1L : 0L;                                           // convert integer rulemod binary to long unsigned integer for use in long logic
    for (ij=0; ij<N2; ij++) {                                               // loop over all sites of 2D torus with side length N
	    i = ij & Nmask;  j = ij >> log2N;                                   // row & column
	    jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
	    ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                         // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;           // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (s=0,k=0,nb1i=0;k<8;k++) {                                      // packs non-zero nb indices in first up to 8*4 bits
            gs=gol[nb[k]];                                                  // whether neighbor is alive
            s += gs;                                                        // s is number of live nbs
            nb1i = (nb1i << (gs<<2)) + (gs*k);                              // nb1i is packed list of live neighbour indices
        }
        if (s>1) {                                                          // if at least 2 neighbours alive
            if (gol[ij]) {                                                  // survival rule in this version just GoL rule
                survive = (s>>2) ? 0L :(s>>1)&1L;                           // 1 if 2 or 3 neighbours are alive : more efficient version of logical (s == 2) || (s == 3)
                newgol[ij]=survive;
                newgolg[ij]=survive*golg[ij];
            }
            else {                                                          // potential birth rule depending on ancestor genes
                if (s==repscheme) {                                         // special rule, repscheme=2 or 4, check next coded nbs & proceed if 3 matching live nbs
                    for (k=kanc=0,s2=0;k<s;k++) {                           // loop only over live neigbours, s2 is number of live nbs in connected 2nd ring
                        nbi = (nb1i>>(k<<2))&0x7;                           // kth live neighbour index
                        ij1 = nb[nbi];                                      // neighbour site ij index
                        gene = golg[ij1];                                   // gene at neighbour site
//                        genematch = gene & 0xffff;                            // 16 bit trailing match subsequence from gene
                        cmask = 0;                                          // connection mask initialized to 0
                        for(l=1;l<4;l++) {                                  // 3 possible connections encoded in 3 16-bit gene words
                            genelink = (gene >> (l<<4)) & 0xffff;           // 16 bit sequences describing possible links
                            if (genelink == 0xffff) cmask = cmask|(1<<(l-1));// set mask only if connection encoded (all ones), later use probs
                        }
                        if(cmask) {                                                  // only if there are some connections
                            i1 = ij1 & Nmask;  j1 = ij1 >> log2N;                                   // row & column
                            j1p1 = ((j1+1) & Nmask)*N; j1m1 = ((j1-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
                            i1p1 =  (i1+1) & Nmask; i1m1 =  (i1-1) & Nmask;                         // toroidal i+1, i-1
                            nb1[0]=j1m1+i1m1; nb1[1]=j1m1+i1; nb1[2]=j1m1+i1p1; nb1[3]=j1*N+i1p1;   //next nbs  0 to 3
                            nb1[4]=j1p1+i1p1; nb1[5]=j1p1+i1; nb1[6]=j1p1+i1m1; nb1[7]=j1*N+i1m1;   //next nbs  4 to 7
                            for(l=1,sl=0;l<4;l++) {
                                if ((cmask>>l)&0x1) {
                                    nbil = ((nbi+6)+l)&0x7; // on the 2nd ring, wrt nb in same direction as k,  the 3 nbs before, at and after (l=1,2,3)
                                    if(gol[nb1[nbil]]) {
//                                        if(genematch == (golg[nb1[nbil]] & 0xffff)) {               // only if next nb genes match
                                            s2 ++;sl ++;
//                                        }
                                    } // if
                                }  // if
                            }  // for
                        } // if
                        if (sl&0x2) kanc = k;   // neighbour contributing 2 or 3 live next shell neighbours serves as ancestor, only works for repscheme==2
                    } // for
                    if(s2==3) {                 // 3 live neighbours in 2nd shell pointed to by live first shell neighbours
                        birth = 1L;
                        newgene = golg[nb[(nb1i>>(kanc<<2))&0x7]];
                    }
                }  // if
                else if (s==3) {                                                             // 3 deterministic choice of ancestor: replication of live neigbour in unique pos
                    for (k=7,nbmask=0L;k>=0;k--) nbmask = (nbmask << 1) + gol[nb[k]];   // compute 8-bit mask of GoL states of 8 neighbours, clockwise starting top left
                    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {                 // compute canonical rotation (minimum) of this mask
                        nbmaskr = ((nbmaskr & 0x1L)<<7) + (nbmaskr>>1);                      // 8 bit rotate right
                        if (nbmaskr < nbmaskrm) {                                          // choose minimal value of mask rotation
                            nbmaskrm = nbmaskr;
                            kmin = k;                                                      // no of times rotated to right
                        }
                    }
                    switch (nbmaskrm) {
                        case 0x07 : k = 1; break;                                  // 00000111
                        case 0x0b : k = 0; break;                                  // 00001011
                        case 0x13 : k = 1; break;                                  // 00010011
                        case 0x19 : k = 0; break;                                  // 00011001
                        case 0x0d : k = 3; break;                                  // 00001101
                        case 0x15 : k = 2; break;                                  // 00010101
                        case 0x25 : k = 3; break;                                  // 00100101
                        default  : {
                            fprintf(stderr,"Error in canonical rotation for three live neighbours \nnbmaskrm = %lx\n",nbmaskrm); k = 0;
                            fprintf(stderr,"Raw Neighbor Pattern: %lx No neighbors %lx\n",
                                    nbmask, gol[nb[0]]+gol[nb[1]]+gol[nb[2]]+gol[nb[3]]+gol[nb[4]]+gol[nb[5]]+gol[nb[6]]+gol[nb[7]]);
                            fprintf(stderr,"\n");
                        }
                    }
                    birth = 1L;
                    newgene = golg[nb[(kmin+k)&0x7]];                               // rotate unique nb k left (kmin) back to orig nb pat
                } // end if (s==3)

                RAND128P(randnr);                                               // inline exp so compiler recognizes auto-vec,
                // compute random events for single bit mutation, as well as mutation position nmut
	            randnr2 = (randnr >> 24) & pmutmask;                                // extract bits from randnr for random trial for 0 on pmutmask
	            r2 = randnr2?0:1;                                                   // 1 if lowest nlog2pmut bits of (bits 24-47 of randnr) are zero, else zero
	            nmut = (randnr >> 48) & 0x3f;                                       // choose mutation position for length 64 gene : from bits 48:53 of randnr
	            // complete calculation of newgol and newgolg, including mutation
	            newgene = newgene ^ (r2*(0x1L<<nmut));                              // introduce single mutation with probability pmut = probmut

	            newgol[ij]  =  gol[ij] | birth ;                                    // new game of life cell value: stays same or set to one from zero if birth
	            newgolg[ij] =  birth*newgene;                                       // if birth (implies empty) then newgene
            }  // end else
        }  // end if s>1
        else {                                                              // else not birth or survival, 0 values
	        newgol[ij]  = 0L;                                                    // new game of life cell value
	        newgolg[ij] = 0L;                                                    // gene dies
        }
        emptysites = emptysites + newgol[ij];                               // keep track of empty sites, same information as total activity of occupied sites
    }  // end for ij

    for (ij=0; ij<N2; ij++) {
	    gol[ij] = newgol[ij];        // copy new gol config to old one
	    golg[ij] = newgolg[ij];      // copy new genes to old genes
    }
}

void genelife_update (int nsteps, int histoflag) {
    /* update GoL for toroidal field which has side length which is a binary power of 2 */
    /* encode without if structures for optimal vector treatment */
    int t;
    long unsigned int *gol, *newgol, *golg, *newgolg;

    static int first = 1;

    if (first) {
        state[0] = rand();state[1] = rand();
        pmutmask = (0x1 << nlog2pmut) - 1;
        first = 0;
    }
    for (t=0; t<nsteps; t++) {
	    gol = planes[curPlane];
	    newgol = planes[newPlane];
	    golg = planesg[curPlane];
	    newgolg = planesg[newPlane];

	    update(gol,golg,newgol,newgolg);
	    if(histoflag) countconfigs();

	    curPlane = (curPlane +1) % numPlane;
	    newPlane = (newPlane +1) % numPlane;
    } /* for t ... */
} /* genelife_update */


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
    int i,j,idx;

    if(N%3 !=0) fprintf(stderr,"Ouch!  Size of offsets array not a multiple of 3 as expected.");
    Noff = N/3;		// Noff global
    if(Noff>24){
        fprintf(stderr,"Too many offsets!  Max=24");
        exit(1);
    }
    numHisto = 1<<Noff;
    histo = (long unsigned int *) calloc(numHisto,sizeof(long unsigned int));

    // install offsets:
    offsets = (int **) calloc(Noff,sizeof(int *));
    for(i=0; i<Noff; i++) offsets[i] = (int *) calloc(3,sizeof(int)); // each containing xoff, yoff, toff.
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

    // initialize plane pointers:
    planes[0] = plane0;planes[1] = plane1;planes[2] = plane2;planes[3] = plane3;
    planes[4] = plane4;planes[5] = plane5;planes[6] = plane6;planes[7] = plane7;
    planesg[0] = planeg0;planesg[1] = planeg1;planesg[2] = planeg2;planesg[3] = planeg3;
    planesg[4] = planeg4;planesg[5] = planeg5;planesg[6] = planeg6;planesg[7] = planeg7;
}

void initialize (int runparams[], int nrunparams, int simparams[], int nsimparams) {
	int ij;

    long unsigned int *gol;
    static unsigned int rmask = (1 << 15) - 1;         // Why 15 bits used here, related to rand(), see next line
    // Range: rand returns numbers in the range of [0, RAND_MAX ), and RAND_MAX is specified with a minimum value of 32,767. i.e. 15 bit

    rulemod = runparams[0];
    repscheme = runparams[1];
    selection = runparams[2];

    nlog2pmut = simparams[0];
    initial1density = simparams[1];
    initialrdensity = simparams[2];
    
    gol = planes[curPlane];
	for (ij=0; ij<N2; ij++) {
		gol[ij] = ((rand() & rmask) < initial1density)?1:0;
	}
}

void initialize_genes (int params[], int nparams) { // params included for possible future use
    int ij,k;
    long unsigned int g;
    long unsigned int *gol;
    long unsigned int *golg;
    static unsigned int rmask = (1 << 15) - 1;

    gol = planes[curPlane];
    golg = planesg[curPlane];

    for (ij=0; ij<N2; ij++) {
        g = 0;
        if (gol[ij] != 0)	{ // if live cell, fill with random genome g or all 1s depending on initialrdensity
            if ((rand() & rmask) < initialrdensity) for (k=0; k<64; k++) g = (g << 1) | (rand() & 0x1);
            else g = 0xf0f0f0f0f0f0f0f0;
        }
        golg[ij] = g;
        if (golg[ij] == 0 && gol[ij] != 0) fprintf(stderr,"zero gene at %d",ij);
    }
    // for (ij=0; ij<40; ij++) fprintf(stderr,"gene at %d %lx\n",ij,golg[ij]);   // test first 40
}

void get_curgol(long unsigned int outgol[], int NN){
    int ij;
    for (ij=0; ij<NN; ij++) {
	    outgol[ij] = planes[curPlane][ij];
    }
}
void get_curgolg(long unsigned int outgolg[], int NN){
    int ij;
    for (ij=0; ij<NN; ij++) {
	    outgolg[ij] = planesg[curPlane][ij];
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
	    fitness = 0;}
        else if (selection == 1) {                                          // non-neutral model with selection for rule departure probability
	    POPCOUNT64C(last, nones);                                     // number of ones in new gene determines fitness
	    fitness = ((nones < 16) ? 0 : (nones - 23));}       // 0 if < 16 otherwise nones-23
        else {                                                              // non-neutral model based on presence of replicase gene
	    POPCOUNT64C(last, nones);                                     // number of ones in new gene determines fitness
	    fitness = ((nones < 16) ? 0 : (nones - 23));}       // 0 if < 16 otherwise nones-23
        fprintf(stderr,"count species %d with gene %lx has counts %lu and %d ones, fitness %d\n",k, golgsc[k][0],golgsc[k][1],nones,fitness);
    }
    printf("cumulative activity = %lu\n",(N2 * (long unsigned int) totsteps) - emptysites);
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

void printxy (long unsigned int gol[],long unsigned int golg[]) {   /* print the game of life configuration */
    int	ij, col, X, Y;
    // https://stackoverflow.com/questions/27159322/rgb-values-of-the-colors-in-the-ansi-extended-colors-index-17-255
    for (ij=0; ij<N2; ij++) {
        if(gol[ij]>0){
	        col = 32+((golg[ij]>>57)&0x7f);
	        X = ij % N;
	        Y = ij / N;
	        printf("%d %d %d ",col,X,Y);
      }
    }
    printf("\n");
}
