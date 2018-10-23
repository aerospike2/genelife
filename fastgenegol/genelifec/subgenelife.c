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
#include <inttypes.h>
#include <time.h>
#include <math.h>

// use Mattias Gustavsson's hashtable (github) for unsigned 64 bit key dictionary
#define HASHTABLE_IMPLEMENTATION
// use 64 bit unsigned key type consistent with this file
#define HASHTABLE_U32 uint32_t
#define HASHTABLE_U64 uint64_t
#define HASHTABLE_SIZE_T uint64_t
#include "hashtable.h"
hashtable_t genetable;
typedef struct genedata {            // value of keys stored for each gene encountered in simulation
            int popcount;            // initialized to 1
            int firstbirthframe;     // initialized to 0
            int lastextinctionframe; // this is initialized to -1, meaning no extinctions yet
            int activity;            // initialized to 0
            int nextinctions;        // initialized to 0
            uint64_t firstancestor;  // this is initialized to a special gene seq not likely ever to occur
            } genedata;
genedata ginitdata = {1,0,-1,0,0,0xfedcba9876543210};  // initialization data structure for gene data
genedata *genedataptr;                      // pointer to a genedata instance
HASHTABLE_SIZE_T const* genotypes;
genedata* geneitems;
#define ASCII_ESC 27                // escape for printing terminal commands, such as cursor repositioning

const int log2N = 9;                // toroidal array of side length N = 2 to the power of log2N
const int N = 0x1 << log2N;
const int N2 = N*N;                 // number of sites in toroidal array
const int Nmask = N - 1;            // bit mask for side length, used instead of modulo operation

int nlog2pmut = 5;                  // pmut (prob of mutation) = 2^(-nlog2pmut), minimum non-zero pmut = 2^-56.
uint64_t pmutmask;                  // binary mask so that prob of choosing zero is pmut, assigned

int totsteps=0;
int rulemod = 1;                    // det: whether to modify GoL rule for 2 and 3 live neighbours : opposite outcome with small probability p0
int selection = 1;                  // fitness of 2 live neighbours: integer value (0), number of ones (1), Norman compete (2), 2 target coding (3)
int repscheme = 1;                  // replication scheme: lowest 3 bits used to define 8 states: bit2 2ndnbs, bit1 not most difft, bit0 2select on 3
int ncoding = 16;                   // maximal distances between sequences for fitness2 == 2
                                    // number of bits used to encode non zero bits in 2nd gene in sequence space for fitness2 == 3
int colorfunction = 0;              // color function choice of hash(0) or functional (1, color classes depends on selection parameter)
int fileinit = 0;                   // 0 input from random field of random or start genes, 1 input from file genepat.dat of 32x32 indexes to start genes
int survival = 0x2;                 // survive mask for two (bit 1) and three (bit 0) live neighbours
int overwritemask = 0x2;            // bit mask for 4 cases of overwrite: bit 0. s==3  bit 1. special birth s==2
int initial1density = (1<<15)>>1;   // initial density of ones in gol as integer value, divide by 2^15 for true density
int initialrdensity = (1<<15)>>1;   // initial density of random genes in live sites, divide by 2^15 for true density
int startgenechoice = 8;            // selection for defined starting genes 0-8 (8 is random 0-7) otherwise choose only particular startgene
uint64_t codingmask;                // mask for ncoding bits
uint64_t  emptysites = 0;           // cumulative number of empty sites during simulation updates

const int startarraysize = 1024;    // starting array size (used when initializing second run)
int arraysize = startarraysize;     // size of trace array (grows dynamically)
int *livesites = NULL;              // dynamic array pointer for statistics of number of live sites over time
int *genestats = NULL;              // dynamic array pointer for statistics of number of 4 genotype classes over time
int *stepstats = NULL;              // dynamic array pointer for statistics of site update types over time

int Noff = 9;                       // number of offsets
int **offsets;                      // array of offsets (2D + time) for planes
int xL=0,xR=0,yU=0,yD=0;            // offsets for border of stats
uint64_t *histo;
int numHisto;

// initialize planes:
#define maxPlane 2
int curPlane = 0;                   // current plane index
int newPlane = 1;                   // new plane index
int numPlane = maxPlane;            // number of planes must be power of 2 to allow efficient modulo plane
uint64_t *planes[maxPlane];         // ring buffer planes of gol array states
uint64_t *planesg[maxPlane];        // ring buffer planes of golg genes
uint64_t plane0[N2];
uint64_t plane1[N2];
uint64_t planeg0[N2];
uint64_t planeg1[N2];
#if maxPlane >= 2
uint64_t plane2[N2];
uint64_t plane3[N2];
uint64_t planeg2[N2];
uint64_t planeg3[N2];
#endif
#if maxPlane >= 4
uint64_t plane4[N2];
uint64_t plane5[N2];
uint64_t plane6[N2];
uint64_t plane7[N2];
uint64_t planeg4[N2];
uint64_t planeg5[N2];
uint64_t planeg6[N2];
uint64_t planeg7[N2];
#endif

uint64_t golgstats[N2];             // 64 bit masks for different events during processing
#define F_notgolrul 0x1
#define F_2_live    0x2
#define F_3_live    0x4
#define F_birth     0x8
#define F_mutation  0x10
#define F_2select   0x20
#define F_survival  0x40
#define F_death     0x80
#define F_golstate  0x100
#define F_3g_same   0x200
#define F_3_livenbs 0xff0000

                                    // Wikipedia "Xorshift" rewritten here as inline macro &
                                    // Vigna, Sebastiano. "xorshift*/xorshift+ generators and the PRNG shootout". Retrieved 2014-10-25.
static uint64_t state[2];           // State for xorshift pseudorandom number generation. The state must be seeded so that it is not zero
#define RAND128P(val) {                                                       \
    uint64_t x = state[0]; uint64_t const y = state[1];                       \
	state[0] = y;	x ^= x << 23;  state[1] = x ^ y ^ (x >> 17) ^ (y >> 26);  \
	val = state[1] + y;}

const uint64_t m1  = 0x5555555555555555; //binary: 0101...           Constants for Hamming distance macro POPCOUNT24C
const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...
#define POPCOUNT64C(x, val) {                /* Wikipedia "Hamming Weight" popcount4c alg */  \
    uint64_t xxxx;                           /* define copy of x argument so that we do not change it */ \
    xxxx = x;                                /* copy x argument */ \
    xxxx -= (xxxx >> 1) & m1;                /* put count of each 2 bits into those 2 bits */ \
    xxxx = (xxxx & m2) + ((xxxx >> 2) & m2); /* put count of each 4 bits into those 4 bits */ \
    xxxx = (xxxx + (xxxx >> 4)) & m4;        /* put count of each 8 bits into those 8 bits */ \
    val = (xxxx * h01) >> 56;}               /* left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... */

void colorgenes(uint64_t gol[],uint64_t golg[], int cgolg[], int NN2) {
    uint64_t gene, gdiff, g2c, mask;
    int ij,d,d2;

    if(colorfunction){
        for (ij=0; ij<NN2; ij++) {
            if (gol[ij]) {
                gene = golg[ij];
                POPCOUNT64C(gene,d);                    // assigns number of ones in gene to d
                switch (selection) {
                        case 0 : mask = ((gene>>40)<<8)+0xff; break;
                        case 1 : mask = ((d+(d<<6)+(d<<12)+(d<<18))<<8) + 0xff; break;
                        case 6 :
                        case 2 : d = d & 0x3; mask = d==3 ? 0xf0f0f0ff : ((0xff<<(d<<3))<<8)+0xff; break;
                        case 3 : g2c = (1L<<ncoding)-1L;gdiff = gene^g2c; POPCOUNT64C(gdiff,d2);
                                 mask = d<d2 ? (d<<26)+0xff : (d2<<10)+0xff; break;
                        case 4 : mask = d < ncoding ? ((0x3f^d)<<20)+0xff : ((64-d < ncoding) ? ((0x3f^d)<<12)+0xff : 0xf0f0f0ff); break;
                        case 5 : mask = d >= 32 ? ((0x3f^(64-d))<<12)+0xff : ((0x3f^d)<<20)+0xff; break;
                        default  : mask = ((d+(d<<6)+(d<<12)+(d<<18))<<8) + 0xff;
                }
                cgolg[ij] = (int) mask;
            }
            else cgolg[ij] = 0;
        }
    }
    else{
        // for(d=0;d<256;d++) counts[d]=0;
        // see https://stackoverflow.com/questions/6943493/hash-table-with-64-bit-values-as-key/33871291
        for (ij=0; ij<NN2; ij++) {
            if (gol[ij]) {
                gene = golg[ij];
                if (gene == 0L) gene = 11778L; // random color for gene==0
                // mask = (gene * 11400714819323198549ul) >> (64 - 8);   // hash with optimal prime multiplicator down to 8 bits
                // mask = (gene * 11400714819323198549ul) >> (64 - 32);  // hash with optimal prime multiplicator down to 32 bits
                mask = gene * 11400714819323198549ul;
                mask = mask >> (64 - 32);   // hash with optimal prime multiplicator down to 32 bits
                mask |= 0x808080ff; // ensure brighter color at risk of improbable redundancy, make alpha opaque
                cgolg[ij] = (int) mask;
            }
            else cgolg[ij] = 0;
        }
    }
    //for(d=0;d<256;d++)
    //    if (counts[d]) fprintf(stderr,"counting hash table hash %d has %d counts\n",d,counts[d]);
}

extern inline void selectone(int s, uint64_t nb2i, int nb[], uint64_t golg[], uint64_t * birth, uint64_t *newgene) {
// birth is returned 1 if ancestors satisfy selection condition. Selection of which of two genes to copy is newgene.
    unsigned int k,d0,d1,d2,d3,dd,swap;                   // number of ones in various gene combinations
    uint64_t livegenes[2],gdiff,gdiff0,gdiff1;                       // various gene combinations
    uint64_t gene2centre;                               // gene function centres in sequence space
    int g0011,g0110,prey,prey2;

    for(k=0;k<2;k++) livegenes[k] = golg[nb[(nb2i>>(k<<2))&0xf]];
    if (selection==0) {                                  // use integer value of sequence as fitness
        *birth = (livegenes[0]^livegenes[1]) ? 1L: 0L;   // birth condition is two genes different
        *newgene = livegenes[0]>livegenes[1] ?  livegenes[0] : livegenes[1]; // choose one with larger gene to replicate
    }
    else {
        POPCOUNT64C(livegenes[0],d0);
        POPCOUNT64C(livegenes[1],d1);
        if (selection==1) {                              // use number of ones in sequence as fitness
            *birth = (d0^d1) ? 1L: 0L;                   // birth condition if two genes different in number of ones
            *newgene= (d0>d1) ? livegenes[0] : livegenes[1];
        }
        else if (selection==2) {                         // use scissors-stone-well-paper game on number ones mod 4
                                                         // scissors 0 stone 1 well 2 paper 3
                                                         // exception to numerical order: sc>pa
            d0=d0&0x3;
            d1=d1&0x3;
            *birth = (d0^d1) ? 1L: 0L;                   // birth if 2 genes differ mod 4 in number of ones
            if(*birth) {
                swap = 0;
                if (d0>d1) {                             // swap d0 and d1 so that smaller one comes first
                    dd = d0;
                    d0 = d1;
                    d1 = dd;
                    swap = 1;
                }
                *newgene = (d0==0 && d1==3) ? livegenes[swap^0] : livegenes[swap^1];
            }
            else *newgene = 0L;
        }
        else if (selection==3) {                         // birth if 2 genes differently functional (Not Yet Working)
            gene2centre = (1L<<ncoding)-1L;              // first ncoding 1s in this sequence
            gdiff  = livegenes[0]^livegenes[1];
            gdiff0 = livegenes[0]^gene2centre;
            gdiff1 = livegenes[1]^gene2centre;
            POPCOUNT64C(gdiff,dd);
            POPCOUNT64C(gdiff0,d2);
            POPCOUNT64C(gdiff1,d3);
            g0011 = d0<dd && d3<dd;
            g0110 = d2<dd && d1<dd;
            *birth = (g0011 != g0110)  ? 1L: 0L;         // birth if 2 genes closer to two different targets than each other
            *newgene= g0011 ? ((d0<d3) ? livegenes[0] : livegenes[1]) : ((d2<d1) ? livegenes[0] : livegenes[1]);
        }
        else if (selection==4) {                         // birth if 2 genes obey 3 distance constraints < ncoding (NYW)
            gdiff=livegenes[0]^livegenes[1];
            POPCOUNT64C(gdiff,dd);
            *birth = (dd<ncoding) && ((d0<ncoding && d1>64-ncoding)|| (d1<ncoding && d0>64-ncoding)) ? 1L: 0L; // birth if 2 genes close enough to targets
            if (d0<ncoding) {if(d0>64-d1) swap=1;else swap=0;}
            else {if(64-d0>d1) swap=1; else swap=0;}
            *newgene= livegenes[swap];
        }
        else if (selection==5) {                         // predator prey model : prey evolves to all 0, predator to all 1
            gdiff=livegenes[0]^livegenes[1];
            gdiff1=livegenes[0]^~(~livegenes[1]);
            POPCOUNT64C(gdiff1,dd);
            prey = d0<32 || d1<32;                       // prey present : newgene is one with less ones, 1 prey : predator wins
            prey2 = d0<32 && d1<32;                       // 2 prey : newgene is one with less ones, 1 prey : predator wins
            *birth = (gdiff && prey && dd<ncoding) ? 1L: 0L;           // birth if different and >=1 prey and close enough match)
            *newgene= (prey2 ? (d0<d1 ? livegenes[0] : livegenes[1]) : (d0<32 ? livegenes[1] : livegenes[0]));
        }
        else if (selection==6) {                         // use next 4 color game on number ones mod 4
                                                         // red 0 green 1 blue 2 white 3
                                                         // exception to numerical order: 0>3 birth only if diff=1
            d0=d0&0x3;
            d1=d1&0x3;
            *birth = ((d0^d1)==1L) ? 1L: 0L;             // birth if 2 genes differ by 1 mod 4 in number of ones
            if(*birth) {
                swap = 0;
                if (d0>d1) {                             // swap d0 and d1 so that smaller one comes first
                    dd = d0;
                    d0 = d1;
                    d1 = dd;
                    swap = 1;
                }
                *newgene = (d0==0 && d1==3) ? livegenes[swap^0] : livegenes[swap^1];
            }
            else *newgene = 0L;
        }
        else fprintf(stderr,"Error: two live gene fitness value %d is not allowed\n",selection);
    }
}

extern inline void selectone_nbs(int s, uint64_t nb2i, int nb[], uint64_t gol[], uint64_t golg[], uint64_t * birth, uint64_t *newgene) {
// birth is returned 1 if ancestors satisfy selection condition. Selection of which of two genes to copy is newgene.
    int k, kanc, l, nb1[8], ij1, i1, j1, j1p1, j1m1, i1p1, i1m1;
    uint64_t nbi, nbil, gene, genelink, s2, sl;
    unsigned int cmask;

    for (k=kanc=0,s2=0;k<s;k++) {                           // loop only over live neigbours, s2 is number of live nbs in connected 2nd ring
        nbi = (nb2i>>(k<<2))&0x7;                           // kth live neighbour index
        ij1 = nb[nbi];                                      // neighbour site ij index
        gene = golg[ij1];                                   // gene at neighbour site
        cmask = 0;                                          // connection mask initialized to 0
        sl = 0;                                             // initialize number of connected live neighbours of this neighbour
        for(l=0;l<2;l++) {                                  // 2 possible connections encoded in 2 16-bit gene words
            genelink = (gene >> ((l+1)<<4)) & codingmask;   // ncoding bit sequences describing possible links: ncoding <=16
            if (genelink == codingmask) cmask = cmask|(1<<l);// set mask only if connection encoded (all ones), later use probs
        }
        if(cmask) {                                                  // only if there are some connections
            i1 = ij1 & Nmask;  j1 = ij1 >> log2N;                                   // row & column
            j1p1 = ((j1+1) & Nmask)*N; j1m1 = ((j1-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
            i1p1 =  (i1+1) & Nmask; i1m1 =  (i1-1) & Nmask;                         // toroidal i+1, i-1
            nb1[0]=j1m1+i1m1; nb1[1]=j1m1+i1; nb1[2]=j1m1+i1p1; nb1[3]=j1*N+i1p1;   //next nbs  0 to 3
            nb1[4]=j1p1+i1p1; nb1[5]=j1p1+i1; nb1[6]=j1p1+i1m1; nb1[7]=j1*N+i1m1;   //next nbs  4 to 7
            for(l=0;l<2;l++) {
                if ((cmask>>l)&0x1) {
                    nbil = (nbi+l)&0x7; // on the 2nd ring, wrt nb in direction k,the 2 nbs at & after (l=0,1)
                    if(gol[nb1[nbil]]) {
                        s2++;sl++;
                    } // if
                }  // if
            }  // for
        } // if
        if (sl&0x2) kanc = k;   // neighbour contributing 2 live next shell neighbours serves as ancestor, only works for repscheme==2
    } // for
    if(s2==3) {                 // 3 live neighbours in 2nd shell pointed to by live first shell neighbours
        *birth = 1L;
        *newgene = golg[nb[(nb2i>>(kanc<<2))&0x7]];
    }
}

void update(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]){
    /* update GoL for toroidal field which has side length which is a binary power of 2 */
    /* encode without if structures for optimal vector treatment */

    int s, k, k1, kmin, nmut;
    int nb[8], nbc, nbch, ij, i, j, jp1, jm1, ip1, im1;
    unsigned int dbirth;
    uint64_t gs, nb1i, nb2i, randnr, randnr2, r2;
    uint64_t nbmask, nbmaskr, nbmaskrm;
    uint64_t newgene, ancestor, livegenes[3];
    uint64_t s2or3, birth, statflag;
    genedata gdata;
    
    totsteps++;
    for (ij=0; ij<N2; ij++) {                                               // loop over all sites of 2D torus with side length N
	    i = ij & Nmask;  j = ij >> log2N;                                   // row & column
	    jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
	    ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                         // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;           // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (s=0L,k=0,nb1i=0;k<8;k++) {                                     // packs non-zero nb indices in first up to 8*4 bits
            gs=gol[nb[k]];                                                  // whether neighbor is alive
            s += gs;                                                        // s is number of live nbs
            nb1i = (nb1i << (gs<<2)) + (gs*k);                              // nb1i is packed list of live neighbour indices
        }
        statflag = 0L;
        s2or3 = (s>>2) ? 0L : (s>>1);                                       // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
        if (s2or3) {                                                        // if 2 or 3 neighbours alive
            if ((s<2)||(s>3)) fprintf(stderr,"s2or3 error s == %d\n",s);
            birth = 0L;
            newgene = 0L;
            if (s&0x1L) {  // s==3                                                 // birth (with possible overwrite)
              statflag |= F_3_live;
              if ((0x1L&overwritemask)|(0x1L&~gol[ij]) ) {                         // central site empty or overwrite mode
                birth = 1L;                                                        // birth flag
                for(k=0;k<s;k++) livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];      // live neighbour genes
                for (k=7,nbmask=0L;k>=0;k--) nbmask = (nbmask << 1) + gol[nb[k]];  // 8-bit mask of GoL states of 8 nbs, clockwise from top left
                statflag |= F_3_livenbs&(nbmask<<16);                              // record live neighbour pattern
                if((livegenes[0]^livegenes[1])|(livegenes[0]^livegenes[2])) {      // genes not all same, need ancestor calculation
                  for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {               // compute canonical rotation (minimum) of this mask
                    nbmaskr = ((nbmaskr & 0x1L)<<7) + (nbmaskr>>1);                // 8 bit rotate right
                    if (nbmaskr < nbmaskrm) {                                      // choose minimal value of mask rotation
                        nbmaskrm = nbmaskr;                                        // neighbor mask rotate min is current rotation
                        kmin = k;                                                  // no of times rotated to right
                    }
                  }
                  if ((repscheme>>1) & 0x1) nbch = nb[kmin];                       // 2,3,6,7... deterministic choice of ancestor: replication of live neigbour in bit 0 of canonical pos
                  else {                                                           // repscheme = 0,1,4,5 for example
                    switch (nbmaskrm) {                //              x07    x0b    x0d    x13    x15    x19    x25
                        case 0x07 : k = 1; break;      // 00000111    |012|  <-
                        case 0x0b : k = 0; break;      // 00001011    |...|  |01.|  <-
                        case 0x0d : k = 3; break;      // 00001101    |...|  |..3|  |0.2|   <-
                        case 0x13 : k = 1; break;      // 00010011           |...|  |..3|  |01.|   <-
                        case 0x15 : k = 2; break;      // 00010101                  |...|  |...|  |0.2|   <-
                        case 0x19 : k = 0; break;      // 00011001                         |..4|  |...|  |0..|   <-
                        case 0x25 : k = 5; break;      // 00100101                                |..4|  |..3|  |0.2|  <-
                        default  : {                   //                                                |..4|  |...|
                                                       //                                                       |.5.|
                            fprintf(stderr,"Error in canonical rotation for three live neighbours \nnbmaskrm = %llx\n",nbmaskrm); k = 0;
                            fprintf(stderr,"Raw Neighbor Pattern: %llx No neighbors %llx\n",
                                nbmask, gol[nb[0]]+gol[nb[1]]+gol[nb[2]]+gol[nb[3]]+gol[nb[4]]+gol[nb[5]]+gol[nb[6]]+gol[nb[7]]);
                            fprintf(stderr,"\n");
                        } //case
                    } //switch
                    nbch = nb[(kmin+k)&0x7];                                // rotate unique nb k left (kmin) back to orig nb pat
                  }
                  if (repscheme & 0x1) {                                    // execute selective replication of one of two otherwise unchosen live genes
                      nb2i = 0L;
                      for(k1=k=0;k<3;k++) {                                 // choice of two other live genes for possible ancestor
                          nbc=(nb1i>>(k<<2))&0x7;
                          if(nb[nbc]!=nbch) nb2i = (nbc<<(k1++<<2))+nb2i;
                      }
                      if ((repscheme>>2)&0x1) selectone_nbs(s,nb2i,nb,gol,golg,&birth,&newgene);
                      else selectone(s,nb2i,nb,golg,&birth,&newgene);
                      if (birth==0L) {                                      // must reset ancestor & birth if no ancestors chosen in selectone
                        dbirth = ((~repscheme)>>3) & 0x1;
                        if(dbirth||rulemod) {
                            newgene = golg[nbch];                           // ALTERED DYNAMICS : less birth for repscheme bit 3 on
                            birth = 1L;
                        }
                      }
                      else statflag |= F_2select;                           // ancestor has been chosen in selectone
                  }
                  else {
                      newgene = golg[nbch];
                  }
                  if (newgene == 0L) fprintf(stderr,"step %d Error with new gene zero: nbmask %llx nbmaskrm %llx kmin %d gol %llx golg %llx newgene %llx ij %d\n",totsteps,nbmask,nbmaskrm,kmin,gol[nb[kmin]],golg[nb[kmin]],newgene,ij);
                } // end if not all live neighbors the same
                else {
                    statflag |= F_3g_same;
                    newgene = livegenes[0];                                // genes all the same : copy first one
                    if((repscheme>>3)&0x1&rulemod) birth = 0L;             // ALTERED DYNAMICS : no birth for 3 identical genes
                }
              }
            }  // end if s==3
            else {  // s==2                                                 // possible birth as exception to GoL rule
                statflag |= F_2_live;
                if (rulemod) {                                              // special rule allowed if rulemod==1
                    if ((0x1L&(overwritemask>>1))|(0x1L&~gol[ij])) {        // either overwrite on for s==2 or central site is empty
                        //for (k=0;k<s;k++) livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]]; // live gene at neighbour site
                        if ((repscheme>>2)&0x1) {
                            selectone_nbs(s,nb1i,nb,gol,golg,&birth,&newgene);
                        }
                        else {    
                            selectone(s,nb1i,nb,golg,&birth,&newgene);
                        }
                        if (birth) statflag |= F_notgolrul|F_2select;
                    }
                }
            }

            if(birth){
                // if (gol[ij]) fprintf(stderr,"birth overwrite event ij %d newgene %llu s %llu\n",ij,newgene,s);
                RAND128P(randnr);                                           // inline exp so compiler recognizes auto-vec,
                // compute random events for single bit mutation, as well as mutation position nmut
                randnr2 = (randnr & pmutmask);                // extract bits from randnr for random trial for 0 on pmutmask
                r2 = ((~pmutmask)||randnr2)?0L:1L;                          // 1 if lowest nlog2pmut bits of randnr are zero, else zero
                nmut = (randnr >> 56) & 0x3f;                               // choose mutation position for length 64 gene : from bits 56:61 of randnr
                // complete calculation of newgol and newgolg, including mutation
                ancestor = newgene;
                newgene = newgene ^ (r2<<nmut);                             // introduce single mutation with probability pmut = probmut
                if(gol[ij]) {                                               // central old gene present: overwritten
                    if((genedataptr = (genedata *) hashtable_find(&genetable, golg[ij])) != NULL) {
                        genedataptr->popcount--;  // if 0 lastextinctionframe updated after whole frame calculated
                    }
                    else fprintf(stderr,"hash storage error 1, gene %llx not stored\n",golg[ij]);
                }

                if((genedataptr = (genedata *) hashtable_find(&genetable, newgene)) != NULL) {
                    genedataptr->popcount++;
                }
                else {
                    gdata=ginitdata;
                    gdata.firstbirthframe = totsteps;
                    gdata.firstancestor = ancestor;
                    hashtable_insert(&genetable, newgene,(genedata *) &gdata);
                }
                newgol[ij]  =  1L;                                          // new game of life cell value: alive
                newgolg[ij] =  newgene;                                     // if birth then newgene
                statflag = statflag | F_birth;
                if (r2) statflag = statflag | F_mutation;
                // if(newgene==0L) fprintf(stderr,"error in writing newgene, previous = %llx, statflag = %llx\n",golg[ij],statflag);
            } // end birth
            else {
                if ((survival&s&0x1L)|((survival>>1)&(~s)&0x1L)|((~rulemod)&0x1L)) {// (surv bit 0 and s==3) or (surv bit 1 and s==2) or not rulemod
                // if ((survival&s&0x1L)|((survival>>1)&(~s)&0x1L)) { // survival bit 0 and s==3, or (survival bit 1 and s==2)
                    newgol[ij]  = gol[ij];                                  // new game of life cell value same as old
                    newgolg[ij] = golg[ij];                                 // gene stays as before, live or not
                    if(gol[ij]) statflag |= F_survival;
                }
                else {
                    if(gol[ij]) {                                           // death : need to update hash table
                        if((genedataptr = (genedata *) hashtable_find(&genetable, golg[ij])) != NULL) {
                            genedataptr->popcount--;
                        }
                        else fprintf(stderr,"hash storage error 2, gene %llx not stored\n",golg[ij]);
                    }
                    newgol[ij]  = 0L;                                       // new game of life cell value dead
                    newgolg[ij] = 0L;                                       // gene dies or stays dead
                    if(gol[ij]) statflag |= F_death|F_notgolrul;
                }
            } // end no birth
        }  // end if s2or3
        else {                                                              // else not birth or survival, 0 values for gol and gene
            if(gol[ij]) {                                                   // death : need to update hash table
                if((genedataptr = (genedata *) hashtable_find(&genetable, golg[ij])) != NULL) {
                    genedataptr->popcount--;
                }
                else fprintf(stderr,"hash storage error 3, gene %llx not stored\n",golg[ij]);
            }
	        newgol[ij]  = 0L;                                               // new game of life cell value
	        newgolg[ij] = 0L;                                               // gene dies
            if(gol[ij]) statflag |= F_death;
        }
        if(gol[ij]) statflag |= F_golstate;
        emptysites = emptysites + newgol[ij];                               // keep track of empty sites, same information as total activity of occupied sites
        golgstats[ij] = statflag;
    }  // end for ij

    for (ij=0; ij<N2; ij++) {       // complete missing hash table records including activities
        if(gol[ij]) {
            if((genedataptr = (genedata *) hashtable_find(&genetable, golg[ij])) != NULL) {
                if(genedataptr->popcount == 0) {
                    genedataptr->lastextinctionframe = totsteps;
                    genedataptr->nextinctions++;
                }
            }
            else fprintf(stderr,"hash storage error 4, gene %llx not stored\n",golg[ij]);
        }
        if(newgol[ij]) {
            if((genedataptr = (genedata *) hashtable_find(&genetable, newgolg[ij])) != NULL)
                genedataptr->activity ++;
            else fprintf(stderr,"hash storage error 5, gene %llx not stored\n",newgolg[ij]);
        }
    }
    
    for (ij=0; ij<N2; ij++) {        // update lattices
	    gol[ij] = newgol[ij];        // copy new gol config to old one
	    golg[ij] = newgolg[ij];      // copy new genes to old genes
    }
}

int statsflag = 1;

void tracestats(uint64_t gol[],uint64_t golg[], uint64_t golgstats[], int NN2) { // trace various stats over time of the simulation
    int ij,cnt,k,d,dc,gt[4],st[10];
    uint64_t gene,statflag;
  
    if (totsteps == arraysize) {                                            // relallocate memory for arrays : double size
        arraysize*=2;
        livesites = (int *)realloc(livesites, arraysize * sizeof(int));
        genestats = (int *)realloc(genestats, arraysize * 4 * sizeof(int));
        stepstats = (int *)realloc(stepstats, arraysize * 10 * sizeof(int));
    }
    for (ij=cnt=0;ij<NN2;ij++) {
        if(gol[ij]) cnt++;
    }
    for(d=0;d<4;d++) gt[d]=0;
    for(k=0;k<10;k++) st[k]=0;
    for (ij=cnt=0;ij<NN2;ij++) {
        statflag = golgstats[ij];
        for(k=0;k<10;k++) if(statflag&(0x1<<k)) st[k]++;
        if(gol[ij]) {
            gene=golg[ij];
            POPCOUNT64C(gene,d);
            switch(selection) {
                case 0: if(d==64) d--; dc=(d>>4)&0x3;break;
                case 1: if(d==64) d--; dc=(d>>4)&0x3;break;
                case 2: dc=d&0x3;break;
                default:if(d==64) d--; dc=(d>>4)&0x3;
            }
            gt[dc]++;
        }
    }

    livesites[totsteps] = cnt;
    for(d=0;d<4;d++) genestats[totsteps*4+d]=gt[d];
    for(k=0;k<10;k++) stepstats[totsteps*10+k]=st[k];
}

void get_stats(int outstats[], int outgtypes[], int outstepstats[], int numStats ){
    int i;
    if(numStats > arraysize){
        fprintf(stderr,"Ack! numStats = %d  > arraysize = %d\n",numStats,arraysize);
        exit(1);
    }
    for(i=0; i<numStats; i++) outstats[i] = livesites[i];
    for(i=0; i<4*numStats; i++) outgtypes[i] = genestats[i];
    for(i=0; i<10*numStats; i++) outstepstats[i] = stepstats[i];
}

void countconfigs(){        // count configs specified by offset array
    // each row of the offset array becomes a bit in an address for the histo array.
    int i,j,k,t,x,y;
    uint64_t *pl, adr, bit;

    for(i=0; i<N; i++){        // rows
        if(i<xL) continue;
        if(i>xR) continue;
        for(j=0; j<N; j++){    // columns
            if(j<yD) continue;
            if(j>yU) continue;
            adr = 0;
            for(k=0; k<Noff; k++){
                x = (j+offsets[k][0]) & Nmask;                          // periodic boundary conditions for N power of 2
                y = (i+offsets[k][1]) & Nmask;
                t = (curPlane+offsets[k][2]) & (numPlane-1);            // periodic in numPlane if this is power of 2
                pl = planes[t];
                bit = *(pl + y*N +x);
                if(bit!= 1L || bit != 0L){                              // check for appropriate unsigned values */
                    fprintf(stderr,"Ack! bit = %llu != 0 or 1\n",bit);
                    exit(1);
                }
                adr = (adr<<1) | bit;
            }
            histo[adr]++;
        }
    }
}

void get_histo(uint64_t outhisto[],int numHistoC){
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

void get_activities(uint64_t actgenes[],int activities[],int ngenesp[]) {
    int k, nlivegenes, nspecies;

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );
    fprintf(stdout,"The number of different species that have ever existed is %d\n",nspecies);
    
    for (k=nlivegenes=0; k<nspecies; k++) {
        if((genedataptr = (genedata *) hashtable_find(&genetable, genotypes[k])) != NULL) {
            if(genedataptr->popcount) {
                actgenes[nlivegenes] = genotypes[k];
                activities[nlivegenes++] = genedataptr->popcount;
            }
        }
        else fprintf(stderr,"get_activities error, no entry for gene %llx in hash table\n", genotypes[k]);
    }
    ngenesp[0] = nlivegenes;
}

void genelife_update (int nsteps, int histoflag) {
    /* update GoL for toroidal field which has side length which is a binary power of 2 */
    /* encode without if structures for optimal vector treatment */
    int t;
    uint64_t *gol, *newgol, *golg, *newgolg;

    for (t=0; t<nsteps; t++) {
        gol = planes[curPlane];
        newgol = planes[newPlane];
        golg = planesg[curPlane];
        newgolg = planesg[newPlane];

        update(gol,golg,newgol,newgolg);
        if(histoflag) countconfigs();
        if(statsflag) tracestats(gol,golg,golgstats,N2);

        curPlane = (curPlane +1) % numPlane;
        newPlane = (newPlane +1) % numPlane;
    } /* for t ... */
} /* genelife_update */

void initialize_planes(int offs[],  int No) {
    int i,j,idx;

    if(No%3 !=0) fprintf(stderr,"Ouch!  Size of offsets array not a multiple of 3 as expected.");
    Noff = No/3;		// Noff global
    if(Noff>24){
        fprintf(stderr,"Too many offsets!  Max=24");
        exit(1);
    }
    numHisto = 1<<Noff;
    histo = (uint64_t *) calloc(numHisto,sizeof(uint64_t));

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
    int tall,tmx = 0;
    int toff, tmn = No;
    for(i=0; i<Noff; i++){
	    toff = offsets[i][2];
	    if(toff>tmx) tmx = toff; // it does not make sense to have +ve values (this would look into future) so tmx = 0
	    if(toff<tmn) tmn = toff;
    }
    if(tmx>0)   {
        fprintf(stderr,"----------------- tmx = %d, tmn = %d",tmx,tmn);
        exit(1);
    }
    tall = tmx-tmn;
    // numPlane = 2 + tall;    // numPlane >= 2
    if(tall>4) numPlane = 8;
    else if (tall > 2) numPlane = 4;
    else numPlane = 2;
    if (numPlane>maxPlane) {
        fprintf(stderr,"Not enough planes defined by maxPlane for given offsets (need > %d)\n",numPlane);
        exit(1);
    }

    // compute xR, xL, yU, yD border offsets from offsets matrix
    int mx;
    int mn;
    int off;
    mn=No; mx=0;
    for(i=0; i<Noff; i++){
	    off = offsets[i][0];	// X
	    if(off>mx) mx = off;
	    if(off<mn) mn = off;
    }
    if(mn<0) xL = -mn; else xL=0;
    if(mx>0) xR = No-mx; else xR=No;
    mn=No; mx=0;
    for(i=0; i<Noff; i++){
	    off = offsets[i][1];	// Y
	    if(off>mx) mx = off;
	    if(off<mn) mn = off;
    }
    if(mn<0) yD = -mn; else yD = 0;
    if(mx>0) yU = No-mx; else yU = No;

    // initialize plane pointers:
    planes[0]  = plane0;  planes[1]  = plane1;
    planesg[0] = planeg0; planesg[1] = planeg1;
#if maxPlane > 2
    planes[2]  = plane2;  planes[3]  = plane3;
    planesg[2] = planeg2; planesg[3] = planeg3;
#endif
#if maxPlane >4
    planes[4]  = plane4;  planes[5]  = plane5;  planes[6]  = plane6;  planes[7]  = plane7;
    planesg[4] = planeg4; planesg[5] = planeg5; planesg[6] = planeg6; planesg[7] = planeg7;
#endif
    curPlane = 0;
    newPlane = 1;
}

int readFile(char * code, char *fileName)
{
  FILE *file;
  char tmp;
  int cnt;
  file = fopen(fileName, "r");
  cnt = -1;
  do
  {
    *code++ = tmp = (char)fgetc(file);
    cnt++;
  } while(tmp != EOF);
  fclose(file);
  return cnt;
}

int writeFile(char *fileName)     // initialize 32x32 genepat file with all empty sites
{
    FILE *file;
    int ij,error;
    file = fopen(fileName, "w");
    error = 0;
    for(ij=0;ij<32*32;ij++) {
        error=fputc(0,file);
        if (error) break;
    }
    fclose(file);
    return error;
}

void initialize (int runparams[], int nrunparams, int simparams[], int nsimparams) {
    int ij,ij1,k,cnt,hcnt,icf;
    uint64_t g;
    uint64_t *gol;
    uint64_t *golg;
    static unsigned int rmask = (1 << 15) - 1;
    static int notfirst = 0;
    // Range: rand returns numbers in the range of [0, RAND_MAX ), and RAND_MAX is specified with a minimum value of 32,767. i.e. 15 bit

    uint64_t startgenes[8];
    
    char *golgin;
    
    srand(1234567);
    state[0] = rand();state[1] = rand();
    cnt = 0;
    
    // writeFile("genepat.dat");

    rulemod = runparams[0];
    repscheme = runparams[1];
    selection = runparams[2];
    overwritemask = runparams[3];
    survival = runparams[4];
    colorfunction = runparams[5];
    fileinit = runparams[6];

    nlog2pmut = simparams[0];
    if(nlog2pmut>56) nlog2pmut=0;                     // need to use top 6-8 bits of 64 bit random nr for position
    pmutmask = (0x1L << nlog2pmut) - 0x1L;            // NB if nlogpmut==0, pmutmask==zero, no mutation.
    initial1density = simparams[1];
    initialrdensity = simparams[2];
    ncoding = simparams[3];
    codingmask = (0x1L<<ncoding)-0x1L;
    startgenechoice = simparams[4];
    
    fprintf(stdout,"runparams %d %d %d %d %d %d %d\n",runparams[0],runparams[1],runparams[2],
                                         runparams[3],runparams[4],runparams[5],runparams[6]);
    fprintf(stdout,"simparams %d %d %d %d %d\n",simparams[0],simparams[1],simparams[2],simparams[3],simparams[4]);
    fprintf(stdout,"pmutmask %llu (NB 0 means no mutation)\n",pmutmask);
    
    switch (selection) {
        case 0: for (k=0;k<4;k++) {startgenes[k]=0xf0f0f0f0f0f0f0f0;startgenes[k+4]=0x0f0f0f0f0f0f0f0f;}; break;
        case 1: for (k=0;k<8;k++) startgenes[k]=((0x1L<<k*3)-1L)<<20;break;
        case 6:
        case 2: for (k=0;k<8;k++) startgenes[k]=(((0x1L<<20)-1L)<<20)+((0x1L<<k)-0x1L);break;
        case 3:
        case 4:
        case 5:
        default: for (k=0;k<8;k++) startgenes[k]=(0x1L<<(4+k*8))-1L;break;
    }

    if ( livesites !=NULL) {
        free(livesites);
        livesites = NULL;
        free(genestats);
        genestats = NULL;
        free(stepstats);
        stepstats = NULL;
    }
    arraysize = startarraysize;
    livesites = (int *) malloc(arraysize * sizeof(int));
    genestats = (int *) malloc(arraysize * 4 * sizeof(int));
    stepstats = (int *) malloc(arraysize * 10 * sizeof(int));

    curPlane = 0;                                           // if we rerun initialize, we want to restart plane cycling from zero
    newPlane = 1;
    gol = planes[curPlane];
    golg = planesg[curPlane];
    
    if(notfirst) hashtable_term(&genetable);
    hashtable_init(&genetable,sizeof(genedata),N2<<2,0);     // initialize dictionary for genes
    notfirst = 1;
    if (fileinit) {           // input from file genepat.dat with max size of 32*32 characters
        golgin = (char *) malloc(32* 32 * sizeof(char));
        icf=readFile(golgin, "genepat.dat");
        if (icf != 32*32) {
            icf = 0;
            fprintf(stderr,"error reading file, %d not 32*32 chars\n",icf);
        }
        for (ij=0; ij<N2; ij++) {
            gol[ij] = 0;
            golg[ij] = 0;
        }
        for (ij1=0; ij1<32*32; ij1++) {
            if(N<32) {fprintf(stderr,"Error, array dim %d too small for file array dim %d\n",N,32);break;}
            ij=(N>>1)-16+(ij1&0x1f)+ N*((N>>1)-16+(ij1>>5));
            if (golgin[ij1] > 0)    {                   // if live cell
                gol[ij] = 1L;
                if(golgin[ij1] <= 8 ) golg[ij] = startgenes[golgin[ij1]-1];
                else if (golgin[ij1]>='0' && golgin[ij1]<'8') golg[ij] = startgenes[golgin[ij1]-'0'];
                else golg[ij] = startgenes[7];
                cnt++;
            }
            else {
                gol[ij] = 0;
                golg[ij] = 0;
            }
            // if (golg[ij] == 0 && gol[ij] != 0) fprintf(stderr,"zero gene at %d\n",ij);
        }

    }
    else {
        for (ij=0; ij<N2; ij++) {
            gol[ij] = ((rand() & rmask) < initial1density)?1:0;
        }
        for (ij=0; ij<N2; ij++) {
            g = 0;
            if (gol[ij] != 0)    { // if live cell, fill with random genome g or randomly chosen startgene depending on initialrdensity
                if ((rand() & rmask) < initialrdensity) for (k=0; k<64; k++) g = (g << 1) | (rand() & 0x1);
                else if (startgenechoice == 8) g = startgenes[rand() & 0x7];
                else if (startgenechoice > 8) fprintf(stderr,"startgenechoice %d out of range\n",startgenechoice);
                else g = startgenes[startgenechoice & 0x7];
                cnt++;
            }
            golg[ij] = g;
            // if (golg[ij] == 0L && gol[ij] != 0L) fprintf(stderr,"zero gene at %d\n",ij);
        }
        // for (ij=0; ij<40; ij++) fprintf(stderr,"gene at %d %llx\n",ij,golg[ij]);   // test first 40
    }

    for (ij=0; ij<N2; ij++) {
        if(gol[ij]) {
            if((genedataptr = (genedata *) hashtable_find(&genetable, golg[ij])) != NULL) {
                genedataptr->popcount++;
            }
            else {
                hashtable_insert(&genetable, golg[ij], (genedata *) &ginitdata);
            }
        }
    }
    // it is possible to enumerate keys and values
    hcnt=hashtable_count(&genetable);
    genotypes = hashtable_keys( &genetable );

    fprintf(stderr,"population size %d with %d different genes\n",cnt,hcnt);
}

void get_curgol(uint64_t outgol[], int NN){
    int ij;
    for (ij=0; ij<NN; ij++) {
	    outgol[ij] = planes[curPlane][ij];
    }
}
void get_curgolg(uint64_t outgolg[], int NN){
    int ij;
    for (ij=0; ij<NN; ij++) {
	    outgolg[ij] = planesg[curPlane][ij];
    }
}

int cmpfunc (const void * pa, const void * pb)
{
   // return ( *(int*)pa - *(int*)pb );
   return ((*(uint64_t*)pa > *(uint64_t*)pb)  ? 1 : -1);
}

int cmpfunc1 ( const void *pa, const void *pb )
{
    const uint64_t *a = pa;
    const uint64_t *b = pb;
    if(a[1] == b[1])
        return a[0] > b[0] ? 1 : -1;
    else
        return (int) (b[1] - a[1]);
}

void countspecies(uint64_t gol[], uint64_t golg[], int N2) {  /* counts numbers of all different species using qsort first */
    int ij, k, ngenes, ijlast, nspecies, counts[N2], nones;
    uint64_t last, golgs[N2], fitness;
    uint64_t golgsc[N2][2];

    for (ij=ngenes=0; ij<N2; ij++) {
        if(gol[ij]) golgs[ngenes++] = golg[ij];                   // only count active sites
        counts[ngenes] = 0;}                                      // initialize sorted gene & count arrays to zero

    qsort(golgs, ngenes, sizeof(uint64_t), cmpfunc);              // sort in increasing gene order
    for (ij=0,k=0,ijlast=0,last=golgs[0]; ij<ngenes; ij++) {      // count each new species in sorted list
        if (golgs[ij] != last) {
            last = golgs[ij];
            counts[k++] = ij - ijlast;
            ijlast = ij;
        }
    }
    counts[k++]=ngenes-ijlast;
    nspecies = k;  // print including 0
    fprintf(stdout,"The number of different species (countspecies) is %d\n",nspecies);

    for (k=0,ij=0;k<nspecies;k++) {     // now condense array to give only different genes with counts
        // printf("species %4d with gene %x has counts %d\n",k, golgs[ij],counts[k]);
        golgs[k]=golgs[ij];
        ij = ij + counts[k];
    }

    for (k=0; k<nspecies; k++) { golgsc[k][0] = golgs[k];  golgsc[k][1] = counts[k];}  // initialize joint gene & count array
    qsort(golgsc, nspecies, sizeof(golgsc[0]), cmpfunc1);                   // sort in decreasing count order


    for (k=1; k<nspecies; k++) {                            // check consistency of hash table data, assuming empty site gene is most frequent
        if((genedataptr = (genedata *) hashtable_find(&genetable, golgsc[k][0])) != NULL) {
                    if(genedataptr->popcount != golgsc[k][1])
                        fprintf(stderr,"popcount %llu <> %d hash error at k = %d\n",golgsc[k][1],genedataptr->popcount,k);
        }
        else fprintf(stderr,"countspecies popcount error, no entry in hash table\n");
    }

    for (k=0; k<nspecies; k++) {
        last = golgsc[k][0];
        POPCOUNT64C(last, nones);
        fitness = 999;
        if (selection == 0) {                                               // 2-live neighbor fitness is integer value
	        fitness = last;
        }
        else if (selection == 1) {                                          // 2-live neighbor fitness is number of ones
	        fitness = (uint64_t) nones;
        }
        else if ((selection == 2)||(selection == 6)) {                      // cyclic 4 species model
	        fitness = nones&0x3;                                            // fitness is species class
        }
        else if (selection == 3){
             fitness = 999;                                                 // undefined, depends on competing sequence
        }
        else if (selection == 4){
             fitness = 999;                                                 // undefined, depends on competing sequence
        }
        else if (selection == 5){
             fitness = 999;                                                 // undefined, depends on competing sequence
        }
        else fprintf(stderr,"selection parameter %d out of range\n",selection);
        fprintf(stdout,"count species %d with gene %llx has counts %llu and %d ones, fitness %llu\n",k, golgsc[k][0],golgsc[k][1],nones,fitness);
    }
    fprintf(stdout,"at step %d cumulative activity = %llu\n",totsteps,(N2 * (uint64_t) totsteps) - emptysites);
    fprintf(stdout,"rulemod\trepscheme\tselection\toverwritemask\tsurvival\n");
    fprintf(stdout,"%d\t%d\t\t%d\t\t%d\t\t%d\n",rulemod,repscheme,selection,overwritemask,survival);
    fprintf(stdout,"nlog2pmut\tinit1\tinitr\tncoding\tstartchoice\n");
    fprintf(stdout,"%d\t\t%d\t%d\t%d\t%d\n",nlog2pmut,initial1density,initialrdensity,ncoding,startgenechoice);
}

int cmpfunc2 (const void * pa, const void * pb)
{
   return ( genotypes[*(int*)pa] > genotypes[*(int*)pb] ? 1 : -1);
}

int cmpfunc3 (const void * pa, const void * pb)
{
   return ( geneitems[*(int*)pa].popcount < geneitems[*(int*)pb].popcount ? 1 : -1);
}

void countspecieshash() {  /* counts numbers of all different species using qsort first */
    int k, golgs[N2], nspecies, nones;
    uint64_t last, fitness;

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );
    fprintf(stdout,"The number of different species is %d\n",nspecies);
    
    for (k=0; k<nspecies; k++) golgs[k] = k;  // initialize sorted genotype array to same order as hash table

    // qsort(golgs, nspecies, sizeof(int), cmpfunc2);                     // sort in increasing gene order
    qsort(golgs, nspecies, sizeof(int), cmpfunc3);                     // sort in decreasing count order

    // for (k=0; k<nspecies; k++) fprintf(stdout,"in countspecieshash genotype %d is %llx\n", k, genotypes[k]);
    for (k=0; k<nspecies; k++) {
        last = genotypes[golgs[k]];
        POPCOUNT64C(last, nones);
        fitness = 999L;
        if (selection == 0) {                                               // 2-live neighbor fitness is integer value
            fitness = last;
        }
        else if (selection == 1) {                                          // 2-live neighbor fitness is number of ones
            fitness = (uint64_t) nones;
        }
        else if ((selection == 2)||(selection == 6)) {                      // cyclic 4 species model
            fitness = nones&0x3;                                            // fitness is species class
        }
        else if (selection == 3){
             fitness = 999L;                                                // undefined, depends on competing sequence
        }
        else if (selection == 4){
             fitness = 999L;                                                // undefined, depends on competing sequence
        }
        else if (selection == 5){
             fitness = 999L;                                                // undefined, depends on competing sequence
        }
        else fprintf(stderr,"selection parameter %d out of range\n",selection);
        if((genedataptr = (genedata *) hashtable_find(&genetable, last)) != NULL) {
            if(genedataptr->popcount)
                fprintf(stdout,"count species %7d with gene %16llx has counts %7d and %3d ones, fitness %llu\n",k,last,
                    genedataptr->popcount,nones,fitness);
        }
        else {
            fprintf(stderr,"countspecieshash popcount error, no entry in hash table\n");
            fprintf(stdout,"count species %d with gene %llx has counts ?? and %d ones, fitness %llu\n",k,last,
                    nones,fitness);
        }
    }
    // fprintf(stdout,"cumulative activity = %llu\n",(N2 * (uint64_t) totsteps) - emptysites);
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

void printxy (uint64_t gol[],uint64_t golg[]) {   /* print the game of life configuration */
    int    ij, col, X, Y;
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
