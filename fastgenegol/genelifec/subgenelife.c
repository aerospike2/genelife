// subgenelife.c
// Written by John S. McCaskill and Norman H. Packard
//
// Project fastgenegol
//
// Created by John McCaskill on 14.07.17.
// Copyright © 2017,2018 European Center for Living Technology. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>
//-----------------------------------------------------------size of array -------------------------------------------------------------------------
const int log2N = 9;                // toroidal array of side length N = 2 to the power of log2N
const int N = 0x1 << log2N;         // only side lengths powers of 2 allowed to enable efficient implementation of periodic boundaries
const int N2 = N*N;                 // number of sites in square-toroidal array
const int Nmask = N - 1;            // bit mask for side length, used instead of modulo operation
const int N2mask = N2 - 1;          // bit mask for array, used instead of modulo operation
//--------------------------------------------------------- main parameters of model ----------------------------------------------------------------
unsigned int rulemod = 1;           // det: whether to modify GoL rule for 2 and 3 live neighbours : opposite outcome with small probability p0
int selection = 1;                  // fitness of 2 live neighbours: 0. integer value   1. number of ones  2. scissors-stone-well-paper: wins over left 1-1-2-2
                                    // 3. scissors-stone-well-paper 1-1-1-1 4. 5. predator prey 6.2 target coding 7. no result
                                    // document 8-13 not yet finished    ******************
unsigned int repscheme = 1;         // replication scheme: lowest 3 bits used to define 8 states: bit2 2ndnbs, bit1 not most difft, bit0 2select on 3
                                    //                     next 2 bits to allow birth failure: bit3 for 3-live-nb and bit4 for 2-live-nb configs
                                    //                     next 2 bits to enableuniquepos case and masking of first neighbours respectively
int survival = 0x2;                 // survive mask for two (bit 1) and three (bit 0) live neighbours
int overwritemask = 0x2;            // bit mask for 4 cases of overwrite: bit 0. s==3  bit 1. special birth s==2
int ncoding = 16;                   // maximal distances between sequences for fitness2 == 2
                                    // number of bits used to encode non zero bits in 2nd gene in sequence space for fitness2 == 3
int nlog2pmut = 5;                  // pmut (prob of mutation) = 2^(-nlog2pmut), minimum non-zero pmut = 2^-56.
uint64_t pmutmask;                  // binary mask so that prob of choosing zero is pmut, assigned
//-----------------------------------------------------------initialization and color parameters------------------------------------------------------
int initial1density = (1<<15)>>1;   // initial density of ones in gol as integer value, divide by 2^15 for true density
int initialrdensity = (1<<15)>>1;   // initial density of random genes in live sites, divide by 2^15 for true density
int startgenechoice = 8;            // selection for defined starting genes 0-8 (8 is random 0-7) otherwise choose only particular startgene
int inittype = 0;                   // 0 input from random field of random or start genes, 1 input from file genepat.dat of 32x32 indexes to start genes
                                    // value i>1: initialized to random field on central ixi block only, outside this zero.
int colorfunction = 0;              // color function choice of 0: hash or 1: functional (color classes depends on selection parameter)
                                    // 2: as in 1 but color sites where last step was non GoL rule yellow, 3: as in 2 but yellow if state produced by non GoL
#define ASCII_ESC 27                /* escape for printing terminal commands, such as cursor repositioning : only used in non-graphic version */
//-----------------------------------------masks for named repscheme bits--------------------------------------------------------------------------------
#define R_0_2sel_3live    0x1       /* 1: employ selection on two least different live neighbours for ancestor with 3 live neighbours */
#define R_1_choose0nb     0x2       /* 1: choose most difft position, 0: choose live neighbour at zero bit in canonical rotation, */
#define R_2_2ndnbs        0x4       /* 1: execute genetically encoded 2nd neighbours of live neighbours to determine birth & ancestor */
#define R_3_enforce3birth 0x8       /* 1: enforce birth for 3 live nbs (with ancestor most difft) in case 2-select fails (req. R_0_2sel_3live==1 for effect) */
#define R_4_enforce2birth 0x10      /* 1: enforce birth for 2 live nbs : choose 1st live nb from top left of canonical config as ancestor (asym!) */
#define R_5_2uniquepos    0x20      /* 1: for 2-live-nb birth, choose clockwise or anticlockwise (for R_1_choose0nb) elt of canonical bunched pair */
#define R_6_checkmasks    0x40      /* 1: check genetically encoded masks to mask out certain live nb pos's: rule as if these not there */
#define R_7_nongolstat    0x80      /* 1: enforce GoL rule if state of central cell was last changed by a non GoL rule */
#define R_8_nongolstatnbs 0x100     /* 1: enforce GoL rule if state of any cell in nbs was last changed by a non GoL rule */
#define R_9_2birth_k0     0x200     /* 1: bit position at start of k1-4 mask for selective subset of 2-births */
#define R_9_12_2birth_k1_4 0x1e00   /* 1: enforce birth for 2 live nbs canonical config for one of k= 1,2,3,4, next 4 bits: choose 1st live nb from TL (asym!) */
//----------------------------------------status flag bits for recording site status in golgstats array---------------------------------------------------
#define F_notgolrul 0x1             /* bit is 1 if last step not a GoL rule*/
#define F_2_live    0x2             /* bit is 1 if exactly 2 live neighbours */
#define F_3_live    0x4             /* bit is 1 if exactly 3 live neighbours */
#define F_birth     0x8             /* bit is 1 if birth (includes overwriting of genes for some parameter values) */
#define F_mutation  0x10            /* bit is 1 if a mutation event occured */
#define F_2select   0x20            /* bit is 1 if the 2 live neighbour selection routine was employed */
#define F_survival  0x40            /* bit is 1 if last step was a 1->1 gol survival */
#define F_death     0x80            /* bit is 1 if last step involved the death of a gene ie 1->0 gol transition */
#define F_golstate  0x100           /* bit is 1 if gol state is 1 */
#define F_golchange 0x200           /* bit is 1 if state changed at last step */
#define F_nongolchg 0x400           /* bit is 1 if state when produced (ie changed to) was made by a non GoL rule */
#define F_3g_same   0x800           /* bit is 1 if exactly 3 live nbs and all 3 have same gene */
#define F_3_livenbs 0xff0000        /* mask for storing configuration of 3 live neighbours : clockwise from top-left neighbour (NW) */
//----------------------------------------------------------hash table implementation of python style dictionary---------------------------------------
#define HASH                        /* commenting out this line removes all hash table stuff */
#ifdef HASH                         // required for countspecieshash() and get_activities(...), otherwise empty routines
                                    // use Mattias Gustavsson's hashtable (github) for unsigned 64 bit key dictionary
#define HASHTABLE_IMPLEMENTATION
                                    // use 64 bit unsigned key type consistent with this file
#define HASHTABLE_U32 uint32_t
#define HASHTABLE_U64 uint64_t
#define HASHTABLE_SIZE_T uint64_t
#include "hashtable.h"              // Gustavsson's file was modified because the 64bit to 32bit key compression code produces 0 for some values
hashtable_t genetable;
typedef struct genedata {           // value of keys stored for each gene encountered in simulation
            int popcount;           // initialized to 1
            int firstbirthframe;    // initialized to 0
            int lastextinctionframe;// this is initialized to -1, meaning no extinctions yet
            int activity;           // initialized to 0
            int nextinctions;       // initialized to 0
            uint64_t firstancestor; // this is initialized to a special gene seq not likely ever to occur for starting genes
            } genedata;
const uint64_t rootgene = 0xfedcba9876543210; // initial special gene as root for genealogies
genedata ginitdata = {1,0,-1,0,0,rootgene};  // initialization data structure for gene data
genedata *genedataptr;              // pointer to a genedata instance
HASHTABLE_SIZE_T const* genotypes;  // pointer to stored hash table keys (which are the genotypes)
genedata* geneitems;                // list of genedata structured items stored in hash table
#endif
//-------------------------------------------------------------------------------------------------------------------------------------------------------
int totsteps=0;                     // total number of simulation steps
int totdisp=0;                      // total number of displayed steps
int statcnts=0;                     // total number of statistic timepts counted
uint64_t codingmask;                // ncoding derived mask for ncoding bits
int nhistG = 0;                     // interval for collecting config histogram data : 0 no collection, nstatG collection with time
int nstatG = 0;                     // interval for collecting other statistical trace data : 0 no collection
int genealogydepth = 0;             // depth of genealogies in current population
//---------------------------------------------------------main arrays in simulation---------------------------------------------------------------------
uint64_t *gol, *golg;               // pointers to gol and golg arrays at one of the plane cycle locations
uint64_t golgstats[N2];             // 64 bit masks for different events during processing
uint64_t newgolgstats[N2];          // 64 bit masks for different events during processing
uint64_t golmix[N2];                // array for calculating coupling of genes between planes for multiplane genelife
uint64_t gene0;                     // uncoupled planes background gene, non zero for selection==10
//------------------------------------------------ arrays for time tracing -----------------------------------------------------------------------------
const int startarraysize = 1024;    // starting array size (used when initializing second run)
int arraysize = startarraysize;     // size of trace array (grows dynamically)
int *livesites = NULL;              // dynamic array pointer for statistics of number of live sites over time
int *genestats = NULL;              // dynamic array pointer for statistics of number of 4 genotype classes over time
int *stepstats = NULL;              // dynamic array pointer for statistics of site update types over time
int *configstats = NULL;            // dynamic array pointer for statistics of gol site configurations (x,y,t) offsets
uint64_t acttrace[N2];              // scrolled trace of last N time points of activity gene trace
int ymax = 1000;                    // activity scale max for plotting : will be adjusted dynamically or by keys
int activitymax;                    // max of activity in genealogical record of current population
uint64_t genealogytrace[N2];        // scrolled trace of genealogies
//------------------------------------------------ planes and configuration offsets----------------------------------------------------------------------
int Noff = 9;                       // number of offsets
int **offsets;                      // array of offsets (2D + time) for planes
int *histo;
int numHisto;
// initialize planes:
#define maxPlane 2                  /* maximum number of planes allowed : values 2,4,8 allowed */
int curPlane = 0;                   // current plane index
int newPlane = 1;                   // new plane index
int numPlane = maxPlane;            // number of planes must be power of 2 to allow efficient modulo plane
uint64_t *planes[maxPlane];         // ring buffer planes of gol array states
uint64_t *planesg[maxPlane];        // ring buffer planes of golg genes
uint64_t plane0[N2];                // gol  0
uint64_t plane1[N2];                // gol  1
uint64_t planeg0[N2];               // golg 0
uint64_t planeg1[N2];               // golg 1
#if maxPlane >= 2
uint64_t plane2[N2];                // gol  2
uint64_t plane3[N2];                // gol  3
uint64_t planeg2[N2];               // golg 2
uint64_t planeg3[N2];               // golg 3
#endif
#if maxPlane >= 4
uint64_t plane4[N2];                // gol  4
uint64_t plane5[N2];                // gol  5
uint64_t plane6[N2];                // gol  6
uint64_t plane7[N2];                // gol  7
uint64_t planeg4[N2];               // golg 4
uint64_t planeg5[N2];               // golg 5
uint64_t planeg6[N2];               // golg 6
uint64_t planeg7[N2];               // golg 7
#endif
//------------------------------------------------------- fast macro random number generator ------------------------------------------------------------
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
//----------------------------------------------------- list of subroutines -----------------------------------------------------------------------------
// colorgenes1          colour genes specified as input parameters
// colorgenes           colour genes specified at current time point
// selectone            select one (or none) of two genes based on selection model parameter selection :  returns birth and newgene
// selectone_nbs        select one of two genes based on pattern of their live 2nd shell neighbours and their genetic encoding
// selectdifft2         select the left or right most difft of two genes bunched with least number of empty genes between them
// selectdifft3         select the unique most different (by symmetry) of three live neighbours
// pack2neighbors       pack all 16 2nd neighbour gol states into a single word bit map, for future use of all 2nd neighbours
// update_gol64         update version for 64 parallel gol planes, currently without genetic coupling (routine not yet used)
// update_gol2          update version for 2 parallel gol planes, coupled by a gene on one plane
// update_gol16         update version for 16 parallel gol planes, coupled by a joint gene for all planes
// update               update the arrays gol, golg, golgstats for a single synchronous time step
// tracestats           record the current stats in time trace
// get_stats            get the traced statistics from python
// countconfigs         count the configs with pthon specified offsets in (x,y,t)
// get_hist             get the histogram from python
// get_activities       get the activity statistics from python
// genelife_update      call update, collect statistics if required and rotate planes
// initialize_planes    initialize periodic sequence of planes to record rolling time window of up to maxPlanes time points (≤8)
// readFile             read file of gol/golg array (32x32) data
// writeFile            write file of gol/golg array (32x32) data
// initialize           initialize simulation parameters and arrays
// get_curgol           get current gol array from python
// get_curgolg          get current golg array from python
// get_acttrace         get current acttrace array from python
// get_curgolgstats     get current golgstats array from python
// cmpfunc              compare gene values as numerical unsigned numbers
// cmpfunc1             compare gene counts in population
// countspecies1        count genes with gene array specified as input parameters
// countspecies         count different genes with genes specified at current time point
// cmpfunc2             compare gene values corresponding to given number index in hash table
// cmpfunc3             compare population counts of hash stored genes
// countspecieshash     count different genes in current population from record of all species that have existed
// activitieshash       calculate array of current activities and update acttrace array of genes in activity plot format
// cmpfunc4             compare birth times of hash stored genes
// genealogies          calculate and display genealogies
// get_sorted_popln_act return sorted population and activities (sorted by current population numbers)
// delay                time delay in ms for graphics
// printxy              terminal screen print of array on xterm
//----------------------------------------------------- begin of subroutines -----------------------------------------------------------------------------
//------------------------------------------------------- colorgenes -------------------------------------------------------------------------------------
void colorgenes1(uint64_t gol[],uint64_t golg[], uint64_t golgstats[], int cgolg[], int NN2) {
    uint64_t gene, gdiff, g2c, mask;
    int ij,d,d2,activity;
    unsigned int color[3],colormax;
    static int numones[16]={0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};

    if(colorfunction==0) { // colorfunction based on multiplicative hash
        // see https://stackoverflow.com/questions/6943493/hash-table-with-64-bit-values-as-key/33871291
        for (ij=0; ij<NN2; ij++) {
            if (gol[ij]) {
                gene = golg[ij];
                if (gene == 0ull) gene = 11778L; // random color for gene==0
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
    else if(colorfunction<4){
        for (ij=0; ij<NN2; ij++) {
            if (gol[ij]) {
                gene = golg[ij];
                POPCOUNT64C(gene,d);                    // assigns number of ones in gene to d
                switch (selection) {
                        case 0 : mask = ((gene>>40)<<8)+0xff; break;
                        case 1 : mask = d==64? 0xffffffff : ((((d&3)<<22)+(((d>>2)&3)<<14)+(((d>>4)&3)<<6))<<8) + 0xff; break;  // number of ones color gradient from black to white
                        case 2 :
                        case 3 : d = d & 0x3; mask = d==3 ? 0xf0f0f0ff : ((0xff<<(d<<3))<<8)+0xff; break;  // scissors-stone-well-paper: red-green-blue-white
                        case 4 : mask = d < ncoding ? ((0x3f^d)<<19)+0xff : ((64-d < ncoding) ? ((0x3f^d)<<11)+0xff : 0xf0f0f0ff); break; // near 0 green, near 1 red, others white
                        case 5 : mask = d >= 32 ? ((0x3f^(64-d))<<11)+0xff : ((0x3f^d)<<19)+0xff; break;  //predators green, prey red
                        case 6 : g2c = (1ull<<ncoding)-1ull;gdiff = gene^g2c; POPCOUNT64C(gdiff,d2);
                                 mask = d<d2 ? (d<<26)+0xff : (d2<<10)+0xff; break;
                        case 7 : g2c = (gene>>8)&((1ull<<ncoding)-1ull);
                                 gdiff = gene&0xff;POPCOUNT64C(gdiff,d2);d = d2>7? 7 : d2;
                                 mask = g2c ? 0xf0f0f0ff : ((0x1f+(d<<5))<<8)+(((gdiff>>4)&0xf)<<27)+(((gdiff&0xf)<<4)<<16)+0xff; break;
                        case 10: for (d=0,mask=0;d<16;d++) {d2=(gene>>(d<<2))&0xf;d2=(d2!=d)?1:0;mask+=(d2*0x2a)<<((d%3)<<3);}
                                mask = (mask<<8)+0x7f7f7fff;break;  // 0x2a is 42, approx 1/6 (ie <16/3) of 255
                        case 11: for (d=0,mask=0;d<16;d++) {d2=(gene>>(d<<2))&0xf;d2=numones[d2];mask+=(d2*0xb)<<((d%3)<<3);}
                                mask = (mask<<8)+0x7f7f7fff;break;  // 0x2a is 42, approx 1/6 (ie <16/3) of 255
                        case 12:
                                d2= (d>>5) ? 5*32 + d-32 : ((d>>4) ? 4*32 + (d-16)*2
                                                         : ((d>>3) ? 3*32 + (d-8)*4
                                                         : ((d>>2) ? 2*32 + (d-4)*8
                                                         : ((d>>1) ? 32 + (d-2)*16
                                                         : d*32))));
                                mask = gene * 11400714819323198549ul;mask = mask >> (64 - 16);
                                mask = ((d2+62)<<8)+(mask<<16)+0xff;break;
                        case 13: mask = ((gene>>40)<<8)+0xff; break;
                        case 14: mask = (d<<20)+((gene&0xff)<<8)+0xff;
                        default  : mask = ((d+(d<<6)+(d<<12)+(d<<18))<<8) + 0xff;
                }
                if(colorfunction==2) {
                    if(golgstats[ij]&F_nongolchg) mask = 0x00ffffff;        // color states changed by non GoL rule yellow
                    if(selection==13) {                                     // color as superposition of biplane gol states
                        if (gol[ij]>>1) {
                            gene = golg[ij];
                            if (gene == 0ull) gene = 11778L; // random color for gene==0
                            mask = gene * 11400714819323198549ul;
                            mask = mask >> (64 - 8);   // hash with optimal prime multiplicator down to 8 bits
                        }
                        else mask = 0;
                        mask = ((0xff*(gol[ij]&1ull))<<24)+((0xff*((gol[ij]>>1)&1ull))<<16)+(mask<<8)+0xff;
                    }
                    else if(selection==12) {                               // color as superposition of multiplane gol states
                        POPCOUNT64C(gol[ij],d);
                        d2= (d>>5) ? 5*32 + d-32 : ((d>>4) ? 4*32 + (d-16)*2
                                                 : ((d>>3) ? 3*32 + (d-8)*4
                                                 : ((d>>2) ? 2*32 + (d-4)*8
                                                 : ((d>>1) ? 32 + (d-2)*16
                                                 : d*32))));
                        mask = gol[ij] * 11400714819323198549ul;mask = mask >> (64 - 16);
                        mask = ((d2+62)<<8)+(mask<<16)+0xff;
                    }
                    else if(selection>=10) {                               // color as superposition of multiplane gol states
                        POPCOUNT64C(gol[ij],d);
                        d2= (d>>3) ? 3*32 + (d-8)*4
                                    : ((d>>2) ? 2*32 + (d-4)*8
                                    : ((d>>1) ? 32 + (d-2)*16
                                    : d*32));
                        mask = gol[ij] * 11400714819323198549ul;mask = mask >> (64 - 16);
                        mask = ((d2+62)<<8)+(mask<<16)+0xff;
                    }
                }
                else if (colorfunction==3) {
                    if(golgstats[ij]&F_notgolrul) mask = 0x00ffffff;  // color states for not GoL rule yellow
                }
                    
                cgolg[ij] = (int) mask;
            }
            else cgolg[ij] = 0;
        }
    }
    else if(colorfunction==4){                    //activities
        for (ij=0; ij<NN2; ij++) {
            gene=acttrace[ij];
            if (gene == rootgene) mask = 0x3f3f3fff;          // grey color for background, all root genes
            else {
                if (gene == 0ull) gene = 11778L; // random color for gene==0
                mask = gene * 11400714819323198549ul;
                mask = mask >> (64 - 32);   // hash with optimal prime multiplicator down to 32 bits
                mask |= 0x808080ff; // ensure brighter color at risk of improbable redundancy, make alpha opaque
            }
            cgolg[ij]= (int) mask;
        }
    }
    else if(colorfunction<8){                    //genealogies
        for (ij=0; ij<NN2; ij++) {
            gene=genealogytrace[ij];
            activity = 0;
            if (gene == rootgene) mask = 0x000000ff;                // black color for root
            else {
                if(colorfunction==7) {
                    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
                        activity = genedataptr->activity;
                    }
                    else {
                        fprintf(stderr,"gene not found in colorfunction for genealogy\n");
                        activity = 0;
                    }
                }
                if (gene == 0ull) gene = 11778L;                    // random color for gene==0
                mask = gene * 11400714819323198549ul;
                mask = mask >> (64 - 32);                           // hash with optimal prime multiplicator down to 32 bits
                if(colorfunction==7) {
                    colormax=0;
                    for(d=0;d<3;d++) if((color[d]=( (mask>>(8+(d<<3))) & 0xff))>colormax) colormax=color[d];
                    for(d=0;d<3;d++) color[d]=color[d]*activity*0xff/(activitymax*colormax);     // rescale colors by activity
                    for(d=0,mask=0xff;d<3;d++) mask |= color[d]<<((d<<3)+8);
                }
                else mask |= 0x808080ff;                                 // ensure brighter color at risk of improbable redundancy, make alpha opaque
            }
            cgolg[ij]=(int) mask;
        }
    }
}

void colorgenes(int cgolg[], int NN2) {
    colorgenes1(gol, golg, golgstats, cgolg, NN2);
}
//------------------------------------------------------- selectone ------------------------------------------------------------
extern inline void selectone(int s, uint64_t nb2i, int nb[], uint64_t golg[], uint64_t * birth, uint64_t *newgene) {
// birth is returned 1 if ancestors satisfy selection condition. Selection of which of two genes to copy is newgene.
    unsigned int k,d0,d1,d2,d3,dd,swap;                  // number of ones in various gene combinations
    uint64_t livegenes[2],gdiff,gdiff0,gdiff1;           // various gene combinations
    uint64_t gene2centre;                                // gene function centres in sequence space
    int g0011,g0110,prey,prey2;

    for(k=0;k<2;k++) livegenes[k] = golg[nb[(nb2i>>(k<<2))&0x7]];
    POPCOUNT64C(livegenes[0],d0);
    POPCOUNT64C(livegenes[1],d1);
    switch (selection) {
        case 0:                                          // integer value of sequence as fitness
            *birth = (livegenes[0]^livegenes[1]) ? 1ull: 0ull; // birth condition is two genes different
            *newgene = livegenes[0]>livegenes[1] ?  livegenes[0] : livegenes[1]; // choose one with larger gene to replicate
            break;

        case 1:                                          // number of ones in sequence as fitness
            *birth = (d0^d1) ? 1ull: 0ull;                   // birth condition is two genes different in number of ones
            *newgene= (d0>d1) ? livegenes[0] : livegenes[1];
            break;
        case 2:                                          // scissors-stone-well-paper game on number ones mod 4
                                                         // scissors 0 stone 1 well 2 paper 3
                                                         // exception to numerical order: sc>pa
            d0=d0&0x3;d1=d1&0x3;
            *birth = (d0^d1) ? 1ull: 0ull;                   // birth if 2 genes differ mod 4 in number of ones
            if(*birth) {
                swap = 0;
                if (d0>d1) {                             // swap d0 and d1 so that smaller one comes first
                    dd = d0; d0 = d1; d1 = dd; swap = 1;
                }
                *newgene = (d0==0 && d1==3) ? livegenes[swap^0] : livegenes[swap^1];
            }
            else *newgene = 0ull;
            break;

        case 3:                                          // 4 color game (next color wins) on number of ones mod 4
                                                         // red 0 green 1 blue 2 white 3
                                                         // exception to numerical order: 0>3 birth only if diff=1
            d0=d0&0x3;
            d1=d1&0x3;
            *birth = ((d0^d1)==1ull) ? 1ull: 0ull;             // birth if 2 genes differ by 1 mod 4 in number of ones
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
            else *newgene = 0ull;
            break;
        case 4:                                          // birth if 2 genes obey 3 distance constraints < ncoding (NYW)
            gdiff=livegenes[0]^livegenes[1];
            POPCOUNT64C(gdiff,dd);
            *birth = (dd<ncoding) && ((d0<ncoding && d1>64-ncoding)|| (d1<ncoding && d0>64-ncoding)) ? 1ull: 0ull; // birth if 2 genes close enough to targets
            if (d0<ncoding) {if(d0>64-d1) swap=1;else swap=0;}
            else {if(64-d0>d1) swap=1; else swap=0;}
            *newgene= livegenes[swap];
            break;
        case 5:                                          // predator prey model: prey-prey evolves to all 0, predator to complement of prey
            gdiff=livegenes[0]^livegenes[1];
            gdiff1=livegenes[0]^(~livegenes[1]);
            POPCOUNT64C(gdiff1,dd);
            prey = (d0<32) || (d1<32);                   // prey present : newgene is one with less ones, 1 prey : predator wins
            prey2 = (d0<32) && (d1<32);                  // 2 prey : newgene is one with less ones, 1 prey : predator wins
            // *birth = (gdiff && prey && dd<ncoding) ? 1ull: 0ull;  // birth if different and >=1 prey and close enough match)
            *birth = ((gdiff && prey2) || (prey && (!prey2) && (dd<ncoding))) ? 1ull: 0ull; // birth if different and >=1 prey and close enough match)
            *newgene= (prey2 ? ((d0<d1) ? livegenes[0] : livegenes[1]) : ((d0<32) ? livegenes[1] : livegenes[0]));
            break;
        case 6:                                         // birth if 2 genes differently functional (Not Yet Working)
            gene2centre = (1ull<<ncoding)-1ull;              // first ncoding 1s in this sequence
            gdiff  = livegenes[0]^livegenes[1];
            gdiff0 = livegenes[0]^gene2centre;
            gdiff1 = livegenes[1]^gene2centre;
            POPCOUNT64C(gdiff,dd);
            POPCOUNT64C(gdiff0,d2);
            POPCOUNT64C(gdiff1,d3);
            g0011 = d0<dd && d3<dd;
            g0110 = d2<dd && d1<dd;
            *birth = (g0011 != g0110)  ? 1ull: 0ull;         // birth if 2 genes closer to two different targets than each other
            *newgene= g0011 ? ((d0<d3) ? livegenes[0] : livegenes[1]) : ((d2<d1) ? livegenes[0] : livegenes[1]);
            break;
        case 7:                                          // no special birth via selection allowed
            *birth = 0ull;
            *newgene = livegenes[0];                     // dummy choice not used
            break;
        default:
            fprintf(stderr,"Error: two live gene fitness value %d is not allowed\n",selection);
            exit(1);
    }
}
//------------------------------------------------------- selectone_nbs -------------------------------------------------------------------------------------
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
        *birth = 1ull;
        *newgene = golg[nb[(nb2i>>(kanc<<2))&0x7]];
    }
}

//------------------------------------------------------- selectdifft2 -------------------------------------------------------------------------------------
extern inline void selectdifft2(uint64_t nbmask, int nb[], uint64_t gol[], unsigned int *kch) {
// selection based on canonical rotation to bunched pair, choose clockwise or anti-clockwise one (for R_1_choose0nb)
    int k,kmin;
    uint64_t nbmaskr, nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {       // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);      // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                            // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                              // neighbor mask rotate min is current rotation
            kmin = k;                                        // no of times rotated to right
        }
    }
    switch (nbmaskrm) {                           //              x03    x05    x09    x11
        case 0x03L : k = 1; break;                // 00000011    |01.|  <-
        case 0x05L : k = 2; break;                // 00000101    |...|  |0.2|  <-
        case 0x09L : k = 3; break;                // 00001001    |...|  |...|  |0..|   <-
        case 0x11ull : k = 4; break;              // 00010001           |...|  |..3|  |0..|   <-
        default  : {                              //                           |...|  |...|
                                                  //                                  |..4|
            fprintf(stderr,"Error in canonical rotation for two live neighbours \nnbmaskrm = %llx\n",nbmaskrm); k = 0;
        } //default case
    } //switch
    if (repscheme & R_1_choose0nb) *kch = kmin;           // replication of live nb in bit 0 of canonical rotation
    else  *kch = (kmin+k)&0x7;                            // replication of live nb in other bit of canonical rotation
}

//------------------------------------------------------- selectdifft3 -------------------------------------------------------------------------------------
extern inline void selectdifft3(uint64_t nbmask, int nb[], uint64_t gol[], unsigned int *kch) {
    int k,kmin;
    uint64_t nbmaskr,nbmaskrm;
    
    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {                // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);             // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                                     // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                       // neighbor mask rotate min is current rotation
            kmin = k;                                                 // no of times rotated to right
        }
    }
    if (repscheme & R_1_choose0nb) *kch = kmin;                       // replication of live neigbour in bit 0 of canonical rotation
    else {                                                            // replication of live neighbour in most different position
        switch (nbmaskrm) {                //              x07    x0b    x0d    x13    x15    x19    x25
            case 0x07L : k = 1; break;      // 00000111    |012|  <-
            case 0x0bL : k = 0; break;      // 00001011    |...|  |01.|  <-
            case 0x0dL : k = 3; break;      // 00001101    |...|  |..3|  |0.2|   <-
            case 0x13L : k = 1; break;      // 00010011           |...|  |..3|  |01.|   <-
            case 0x15L : k = 2; break;      // 00010101                  |...|  |...|  |0.2|   <-
            case 0x19L : k = 0; break;      // 00011001                         |..4|  |...|  |0..|   <-
            case 0x25L : k = 5; break;      // 00100101                                |..4|  |..3|  |0.2|  <-
            default  : {                   //                                                |..4|  |...|
                                           //                                                       |.5.|
                fprintf(stderr,"Error in canonical rotation for three live neighbours \nnbmaskrm = %llx\n",nbmaskrm); k = 0;
                fprintf(stderr,"Raw Neighbor Pattern: %llx No neighbors %llx\n",
                                nbmask, gol[nb[0]]+gol[nb[1]]+gol[nb[2]]+gol[nb[3]]+gol[nb[4]]+gol[nb[5]]+gol[nb[6]]+gol[nb[7]]);
                fprintf(stderr,"\n");
            } //default case
        } //switch
        *kch = (kmin+k)&0x7;                                        // rotate unique nb k left (kmin) back to orig nb pat
    }
}

//------------------------------------------------------- hash gene inline fns ----------------------------------------------------------------------
extern inline void hashaddgene(uint64_t gene,uint64_t ancestor) {
#ifdef HASH
        genedata gdata;
        if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) genedataptr->popcount++;
        else {
            gdata=ginitdata;
            gdata.firstbirthframe = totsteps;
            gdata.firstancestor = ancestor;
            hashtable_insert(&genetable, gene,(genedata *) &gdata);
        }
#endif
}

extern inline void hashdeletegene(uint64_t gene,char errorformat[]) {
#ifdef HASH
        if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) genedataptr->popcount--;
        else fprintf(stderr,errorformat,totsteps,gene);
#endif
}

extern inline void hashreplacegene(uint64_t gene1,uint64_t gene2,uint64_t ancestor,char errorformat[]) {
#ifdef HASH
        genedata gdata;
        if((genedataptr = (genedata *) hashtable_find(&genetable, gene1)) != NULL) genedataptr->popcount--;
        else fprintf(stderr,errorformat,totsteps,gene1);
        if((genedataptr = (genedata *) hashtable_find(&genetable, gene2)) != NULL) genedataptr->popcount++;
        else {
            gdata=ginitdata;
            gdata.firstbirthframe = totsteps;
            gdata.firstancestor = ancestor;
            hashtable_insert(&genetable, gene2,(genedata *) &gdata);
        }
#endif
}

extern inline void hashgeneextinction(uint64_t gene,char errorformat[]) {
#ifdef HASH
        if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
            if(genedataptr->popcount == 0) {
                genedataptr->lastextinctionframe = totsteps;
                genedataptr->nextinctions++;
            }
        }
        else fprintf(stderr,errorformat,3,totsteps,gene);
#endif
}

extern inline void hashgeneactivity(uint64_t gene,char errorformat[]) {
#ifdef HASH
        if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) genedataptr->activity ++;
        else fprintf(stderr,errorformat,4,totsteps,gene);
#endif
}

//------------------------------------------------------- pack2neighbors ----------------------------------------------------------------------------------
void pack2neighbors(uint64_t gol[], uint64_t golp[]) {              // routine for future use of all 2nd neighbours
#define deltaxy(ij,x,y)  (ij - (ij&Nmask) + ((ij+(x)&Nmask) + (y)*N)) & N2mask
    unsigned int ij,k,s;
    uint64_t gs;
    //int nb1x[8] = {-1,0,1,1,1,0,-1,-1};
    //int nb1y[8] = {-1,-1,-1,0,1,1,1,0};
    int nb2x[16] = {-2,-1,0,1,2,2,2,2,2,1,0,-1,-2,-2,-2,-2};
    int nb2y[16] = {-2,-2,-2,-2,-2,-1,0,1,2,2,2,2,2,1,0,-1};

    for (ij=0;ij<N2;ij++) golp[ij] = 0ull;
    for (ij=0;ij<N2;ij++) {
        for(k=s=0;k<16;k++) {
            gs=gol[deltaxy(ij,nb2x[k],nb2y[k])];
            // golp[ij] |= gs<<(k+1);
            golmix[ij] |= (k*gs) << (s<<2);
            s+=gs;
        }
    }
}

//------------------------------------------------------- update_gol64 -----------------------------------------------------------------------------------
void update_gol64(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]) {  // routine for 64x-gol packed update, no plane coupling
#define deltaxy(ij,x,y)  (ij - (ij&Nmask) + ((ij+(x)&Nmask) + (y)*N)) & N2mask
    unsigned int ij,k;
    uint64_t gs,newgs,sums00,sums10,sums01,sums11,sums02,sums12,sums16[4];
    uint64_t s,su,sg3,s3,s2or3,sgol;
    int nb1x[8] = {-1,0,1,1,1,0,-1,-1};
    int nb1y[8] = {-1,-1,-1,0,1,1,1,0};
    const uint64_t r1= 0x1111111111111111;
    const uint64_t r3= 0x3333333333333333;
    const uint64_t r5= 0x5555555555555555;
    const uint64_t ra= 0xaaaaaaaaaaaaaaaa;
    const uint64_t rc= 0xcccccccccccccccc;
    
    totsteps++;
    if(!(totsteps%10)) fprintf(stderr,"iteration step %d\r",totsteps);
    for (ij=0; ij<N2; ij++) {                                               // loop over all sites of 2D torus with side length N
        sums00 = sums10 = sums01 = sums11 = sums02 = sums12 = 0ull;
        for(k=0;k<3;k++) {
            gs=gol[deltaxy(ij,nb1x[k],nb1y[k])];
            sums00 += gs&r5;
            sums10 += (gs&ra)>>1;
        }
        for(k=3;k<6;k++) {
            gs=gol[deltaxy(ij,nb1x[k],nb1y[k])];
            sums01 += gs&r5;
            sums11 += (gs&ra)>>1;
        }
        for(k=6;k<8;k++) {
            gs=gol[deltaxy(ij,nb1x[k],nb1y[k])];
            sums02 += gs&r5;
            sums12 += (gs&ra)>>1;
        }
        sums16[0] = (sums00&r3)+(sums01&r3)+(sums02&r3);
        sums16[1] = (sums10&r3)+(sums11&r3)+(sums12&r3);
        sums16[2] = ((sums00&rc)>>2)+((sums01&rc)>>2)+((sums02&rc)>>2);
        sums16[3] = ((sums10&rc)>>2)+((sums11&rc)>>2)+((sums12&rc)>>2);
        newgs=0;
        for (k=0;k<4;k++) {
            s = sums16[k];su = s&rc;
            sg3=(((su>>1)|su)>>2) & r1;
            s2or3 = (~sg3)&(s>>1)&r1;
            s3=s2or3&s&r1;
            sgol=(s2or3&(gol[ij]|s3))&r1;
            newgs|=sgol<<k;
        }
        //for (k=0;k<64;k++) {                                                    // explicit loop over 64 bits
        //    s = (sums16[k&0x3]>>(k>>2))&0xf;
        //    s2or3 = (s>>2) ? 0ull : (s>>1);                                     // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
        //    newgs |= (s2or3 & (gol[ij] |(s&0x1ull)))<<k;                        // GoL standard calculation next state
        //}
        newgol[ij]=newgs;
        newgolg[ij]=golg[ij];                                                     // currently no gene changes
    }
    for (ij=0; ij<N2; ij++) {        // update lattices
        gol[ij] = newgol[ij];        // copy new gol config to old one
        golg[ij] = newgolg[ij];      // copy new genes to old genes
    }
}

//------------------------------------------------------- update_gol2 -----------------------------------------------------------------------------------
void update_gol2(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]) {
// gol states on one plane interact with genes hard coupled to the gol states on a second plane
// this may also be viewed simply as genes with gol dynamics interacting with non-genetic entities with gol dynamics
// both gol planes are coded as two bottom bits in the gol[i,j] array
// coupling between the two planes occurs through complementary (1<->0) exact matching of the local Moore neighborhoods (8-sites)
// survival and/or birth rules are negated (opposite outcome) when coupling is established
#define deltaxy(ij,x,y)  (ij - (ij&Nmask) + ((ij+(x)&Nmask) + (y)*N)) & N2mask
    unsigned int ij,ij1,k,kch,nmut;
    uint64_t s0,s1,s0_2or3,s1_2or3,nbmask0,nbmask1,ancestor,newgene;
    uint64_t golmix,gol0,gol1,newgol1;
    uint64_t randnr, randnr2, rand2, statflag;
    int nb[8];
    int nb1x[8] = {-1,0,1,1,1,0,-1,-1};
    int nb1y[8] = {-1,-1,-1,0,1,1,1,0};


    totsteps++;
    if(!(totsteps%10)) fprintf(stderr,"iteration step %d\r",totsteps);
    
    for (ij=0; ij<N2; ij++) {                                                   // loop over all sites of 2D torus with side length N
        for(k=0,s0=s1=nbmask0=nbmask1=0ull;k<8;k++) {                           // compute no of live neighbours using golmix (coupled plane neighbours)
            ij1=nb[k]=deltaxy(ij,nb1x[k],nb1y[k]);
            gol0=gol[ij1]&1ull;
            gol1=(gol[ij1]>>1)&1ull;
            s0 += gol0;
            s1 += gol1;
            nbmask0 |= (gol0<<k);
            nbmask1 |= (gol1<<k);
        }
        gol0=gol[ij]&1ull;
        gol1=(gol[ij]>>1)&1ull;
        if(rulemod && (nbmask0==nbmask1)) golmix=1ull;                          // also tried ==~nbmask1 but virtually no coupling
        else golmix=0ull;
        s0_2or3 = (s0>>2) ? 0ull : (s0>>1);                                     // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
        s1_2or3 = (s1>>2) ? 0ull : (s1>>1);                                     // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
        newgol[ij] = s0_2or3 ? (gol0 ? 1ull : (s0&1ull ? 1ull : golmix )) : 0ull;  // GoL calculation next state for non-genetic gol plane
        newgol1    = s1_2or3 ? (gol1 ? 1ull : (s1&1ull ? 1ull : golmix )) : 0ull;  // GoL calculation next state for genetic gol plane
        newgol[ij] |= (newgol1<<1);                                               // combine gol planes
        statflag = 0ull;
        if(golmix&(s0_2or3|s1_2or3)) statflag |= F_notgolrul;                   // not gol rule
        if (newgol1&(~gol1)&1ull) {                                             // birth in gene plane without overwrite
            if(s1&1ull)    // 3 nbs
                selectdifft3(nbmask1, nb, gol, &kch);
            else      // 2 nbs
                selectdifft2(nbmask1, nb, gol, &kch);
            ij1=deltaxy(ij,nb1x[kch],nb1y[kch]);
            newgene = golg[ij1];
            RAND128P(randnr);                                                   // inline exp so compiler recognizes auto-vec,
            randnr2 = (randnr & pmutmask);                                      // extract bits from randnr for random trial for 0 on pmutmask
            rand2 = (!pmutmask||randnr2)?0ull:1ull;                             // 1 if lowest nlog2pmut bits of randnr are zero, else zero
            nmut = (randnr >> 56) & 0x3f;                                       // choose mutation pos for length 64 gene : from bits 56:61 of randnr
            ancestor = newgene;
            newgene = newgene ^ (rand2<<nmut);                                  // introduce single mutation with probability pmut = probmut
            newgolg[ij]=newgene;
            statflag = statflag | F_birth;
            if (rand2) statflag = statflag | F_mutation;
            if(statflag&F_notgolrul) statflag |= F_nongolchg;
            hashaddgene(newgene,ancestor);
        } //end if birth
        else if (newgol1) {                                                     // gene survival
            newgolg[ij]=golg[ij];
            statflag = statflag | F_survival;
        }
        else {                                                                  // gene death if alive
            if(gol1) {                                                          // site with gene
                hashdeletegene(golg[ij],"step %d hash storage error 2 in update_gol2, gene %llx not stored\n");
                newgolg[ij]=gene0;
                statflag = statflag | F_death;
            }
            else newgolg[ij]=golg[ij];                                          // stays dead
        } // end else
        newgolgstats[ij] = statflag;
    } // for ij

    for (ij=0; ij<N2; ij++) {       // complete missing hash table records including activities
        if(gol[ij]>>1) hashgeneextinction(golg[ij],"hash storage error 3 in update_gol2, gene %llx not stored\n");
        if(newgol[ij]>>1) hashgeneactivity(newgolg[ij],"hash storage error 4 in update_gol2, gene %llx not stored\n");
    }

    for (ij=0; ij<N2; ij++) {        // update lattices
        gol[ij] = newgol[ij];        // copy new gol config to old one
        golg[ij] = newgolg[ij];      // copy new genes to old genes
        golgstats[ij] = newgolgstats[ij]; // copy new golg stat flags to old genes
    }
}

//------------------------------------------------------- update_gol16 -----------------------------------------------------------------------------------
void update_gol16(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]) {
// genes specify on which of 16 planes the gol states at that site are visible: specified by 16 bits 0,4,8,... in gene
// when multiple copy processes are active for gol states, the one with the lowest bit position is applied to most difft gene copying
// gol states looked up for neighbors are the OR of the current plane gol value with gol values on planes specified by the gene
// optionally genes also specify a unique copy plane for a gene: it can only be copied by a rule active on this plane
// in this optional case, copy plane could be lowest non zero bit in gene at pos m for which m%4 == 1 (if none then gene is not copied)
#define deltaxy(ij,x,y)  (ij - (ij&Nmask) + ((ij+(x)&Nmask) + (y)*N)) & N2mask
    unsigned int ij,ij1,k,kmin,p,p1,nmut,debcnt,mask,np,npmask;
    uint64_t s,sm,sm3,s3sm3,su,s3,sg3,s2or3,nbmask,nbmaskr,nbmaskrm,ancestor,newgene,golsh,pmask;
    uint64_t *golm;
    //uint64_t s0,sm;
    uint64_t randnr, randnr2, rand2, statflag;
    int nb1x[8] = {-1,0,1,1,1,0,-1,-1};
    int nb1y[8] = {-1,-1,-1,0,1,1,1,0};
    int npmasks[17] = {1,1,3,3,7,7,7,7,0xf,0xf,0xf,0xf,0xf,0xf,0xf,0xf,0xf};
    const uint64_t r1= 0x1111111111111111;
    const uint64_t rc= 0xcccccccccccccccc;

    totsteps++;
    
    np = (repscheme && repscheme<=16) ? repscheme : 16;                      // number of planes used, 16 if 0
    npmask = npmasks[np-1];                                                 // mask for other plane index
    pmask = np==16 ? 0xffffffffffffffff : (1ull<<(np<<2))-1;                // mask for bits included in np planes
    
    if(!(totsteps%10)) fprintf(stderr,"iteration step %d\r",totsteps);
    if(rulemod) {
      if (selection == 10) {                                                 // gene codes for one secondary plane for each plane  (OLD 10)
        for (ij=0; ij<N2; ij++) {                                            // loop over all sites of 2D torus with side length N
            for (p=0,golsh=0ull;p<np;p++) {                                  // 16 different planes, one every 4 bits in gol,golg
                p1 = (golg[ij]>>(p<<2))&npmask;
                p1 &= ncoding ? npmasks[p+1] : 0xf;                          // 1st implementation of Norman's idea to grow plane community
                golsh |= ((gol[ij]>>(p1<<2))&1ull)<<(p<<2);                  // collect active bits from 2nd planes pointed to by gene
            }
            golsh &= pmask;
            if (survival&0x1)   golmix[ij] = golsh|gol[ij];                  // bits active from secondary plane or primary plane
            else golmix[ij] = golsh;                                         // bits active from secondary plane only
        }
      }
      else { // selection == 11                                              // gene codes for 4(2) nearest plane neighbour masks
        for (ij=0; ij<N2; ij++) {                                            // loop over all sites of 2D torus with side length N
            for (p=0,golsh=0ull;p<np;p++) {                                  // 16 different planes, one every 4 bits in gol,golg
                mask = (golg[ij]>>(p<<2))&0xf;                               // bit mask indicating which nearby planes visible from plane
                for (k=1;k<3;k++) {                                          // 2 nearest planes variant
                //for (k=0;k<4;k++) {                                        // 4 nearest planes variant
                    if((mask>>k)&0x1) {
                        p1 = (p+k+(k>>1)-2)&npmask;
                        golsh ^= ((gol[ij]>>(p1<<2))&0x1ull)<<(p<<2);        // nearby planes p1=p+{-2,-1,1,2} for k 0,1,2,3   or p+{-1,1} for k 1,2
                    }
                }
            }
            golsh &= pmask;
            if (survival&0x1)   golmix[ij] = golsh|gol[ij];                  // bits active from secondary plane or primary plane
            else golmix[ij] = golsh;                                         // bits active from secondary plane only
        }
      }
    }
    else { // no gol rule modification, independent planes
         for (ij=0; ij<N2; ij++) golmix[ij] = gol[ij];                       // loop over all sites of 2D torus with side length N
    }
    // for (ij=0; ij<N2; ij++) if ((golmix[ij]&r1) != golmix[ij]) fprintf(stderr,"error in golmix calculation %llx",golmix[ij]);
    
    for (ij=0; ij<N2; ij++) {                                                   // loop over all sites of 2D torus with side length N
        for(k=0,s=0ull,sm=0ull;k<8;k++) {                                       // compute no of live neighbours using golmix (coupled plane neighbours)
            ij1=deltaxy(ij,nb1x[k],nb1y[k]);
            sm +=golmix[ij1];
            s +=gol[ij1];
        }
        sm3 = (~sm>>3)&(~sm>>2)&(sm>>1)&sm&r1;                                  // sum of golmix bits is exactly 3, calculated for each plane
        su = s&rc;                                                              // upper 2 bits (3,2) of sum of live neighbours : non zero if sum 4-8
        sg3 = (((su>>1)|su)>>2) & r1;                                           // for each plane 1 if sum is gt 3
        s2or3 = (~sg3)&(s>>1)&r1;                                               // sum is 2 or 3 for each plane, ie not gt 3 and bit 2 is 1
        s3 = s2or3&s;                                                           // sum is 3 (s2or3 and bit 0 is 1) calculated for each plane
        s3sm3 = s3|sm3;
        if (survival&0x2)
            newgol[ij]=s2or3&(gol[ij]|s3sm3);                                  // parallel gol rule with coupled sum s for each plane
        else
            newgol[ij]=s2or3&(s3sm3);                                          // no survival only birth
        statflag = 0ull;
        //if((golmix[ij]^gol[ij])&s2or3) statflag |= F_notgolrul;               // not gol rule in at least one plane
        if((s3^sm3)&s2or3) statflag |= F_notgolrul;                             // not gol rule in at least one plane
        if (s3sm3) {                                                           // birth in at least one plane, implement for lowest such birth plane
          for(p=0;p<np;p++) {
            if((s3sm3>>(p<<2))&1ull) {
                if((s3>>(p<<2))&1ull) golm = gol;
                else golm = golmix;
                for(k=0,nbmask=0ull;k<8;k++) {
                    ij1=deltaxy(ij,nb1x[k],nb1y[k]);
                    nbmask |= ((golm[ij1]>>(p<<2))&1ull)<<k;
                }
                POPCOUNT64C(nbmask,debcnt);
                if(debcnt!=3) fprintf(stderr,"Error with nbmask %llx at p %d with s3 %llx s %llx\n",nbmask,p,s3,s);
                for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {               // compute canonical rotation (minimum) of this mask
                    nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);            // 8 bit rotate right
                    if (nbmaskr < nbmaskrm) {                                    // choose minimal value of mask rotation
                        nbmaskrm = nbmaskr;                                      // neighbor mask rotate min is current rotation
                        kmin = k;                                                // no of times rotated to right
                    }
                }
                switch (nbmaskrm) {
                    case 0x07L : k = 1; break;                                   // 00000111
                    case 0x0bL : k = 0; break;                                   // 00001011
                    case 0x0dL : k = 3; break;                                   // 00001101
                    case 0x13L : k = 1; break;                                   // 00010011
                    case 0x15L : k = 2; break;                                   // 00010101
                    case 0x19L : k = 0; break;                                   // 00011001
                    case 0x25L : k = 5; break;                                   // 00100101
                    default  : fprintf(stderr,"Error in canonical rotation for three live neighbours \nnbmaskrm = %llx\n",nbmaskrm); k = 0;
                }
                k=(kmin+k)&0x7;
                ij1=deltaxy(ij,nb1x[k],nb1y[k]);
                newgene = golg[ij1];
                RAND128P(randnr);                                                 // inline exp so compiler recognizes auto-vec,
                randnr2 = (randnr & pmutmask);                                    // extract bits from randnr for random trial for 0 on pmutmask
                rand2 = (!pmutmask||randnr2)?0ull:1ull;                           // 1 if lowest nlog2pmut bits of randnr are zero, else zero
                nmut = (randnr >> 56) & (0x3|(npmask<<2));                        // choose mutation pos for length np*4 gene : from bits 56:61 of randnr
                ancestor = newgene;
                newgene = newgene ^ (rand2<<nmut);                                // introduce single mutation with probability pmut = probmut
                newgolg[ij]=newgene;
                statflag = statflag | F_birth;
                if (rand2) statflag = statflag | F_mutation;
                newgolgstats[ij] = statflag;
                hashreplacegene(golg[ij],newgene,ancestor,"step %d hash storage error 1 in update_gol16, gene %llx not stored\n");
                p=np;                                                             // break from loop at first birth
            } //if 3 live nbs at p
        } // for p
      }  // if (s3)
      else {
        newgolg[ij]=gene0;
        if(golg[ij]!= gene0) {
            newgene= gene0;                                                     // default is relax to uncoupled gene
            hashreplacegene(golg[ij],newgene,rootgene,"step %d hash storage error 1 in update_gol16, gene %llx not stored\n");
        }
      } // end else (s3)
    } // for ij
    
    for (ij=0; ij<N2; ij++) {       // complete missing hash table records including activities, NB all sites have genes in gol16 but don't record activity for zerogene
        hashgeneextinction(golg[ij],"hash storage error 4 in update_gol16, gene %llx not stored\n");
        if(newgolg[ij]!=gene0) hashgeneactivity(newgolg[ij],"hash storage error 5 in update_gol16, gene %llx not stored\n");
    }

    for (ij=0; ij<N2; ij++) {        // update lattices
        gol[ij] = newgol[ij];        // copy new gol config to old one
        golg[ij] = newgolg[ij];      // copy new genes to old genes
        golgstats[ij] = newgolgstats[ij]; // copy new golg stat flags to old genes
    }


}
//------------------------------------------------------- update_lut_sum -----------------------------------------------------------------------------------
void update_lut_sum(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]){
    /* update GoL for toroidal field which has side length which is a binary power of 2 */
    /* encode without if structures for optimal vector treatment */
    /**********
        0 <= s <= 8
        genome =
        8 bits for each possible s value for birth (center site 0)
        9 bits for each possible s value for survival/death (center site 1)
        (assume s = 0 => 0 always)
        using b0 for s=0, b1 for s=1, etc
        GOL = 00001000|000001100
            = 0 0001 0000 0000 1100
            = 0x100c
    ***********/

    int s, k, nmut, idx;
    unsigned int survivemask;
    int nb[8],  ij, i, j, jp1, jm1, ip1, im1;
    uint64_t gene, gols, nb1i, randnr, randnr2, r2;
    uint64_t newgene, ancestor, livegenes[8];
    uint64_t  birth, statflag, ncodingmask;
    
    totsteps++;
    survivemask=(repscheme>>16)&0x1ff;                                      // 9 bits of repscheme 16-24 used to limit space of rules for survival
    if (ncoding<1 || ncoding>7) ncoding = 7;                                // ncoding out of range gives maximum currently allowing value of ncoding
    ncodingmask = (1ull<ncoding)-1ull;                                      // mask for number of bits coding for each lut rule: <=7 for survival, <=3 for b and s

    for (ij=0; ij<N2; ij++) {                                               // loop over all sites of 2D torus with side length N
        gene = golg[ij];
        i = ij & Nmask;  j = ij >> log2N;                                   // row & column
        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                         // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;           // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (s=0,k=0,nb1i=0;k<8;k++) {                                     // packs non-zero nb indices in first up to 8*4 bits
            gols=gol[nb[k]];                                                // whether neighbor is alive
            s += gols;                                                      // s is number of live nbs
            nb1i = (nb1i << (gols<<2)) + (gols*k);                          // nb1i is packed list of live neighbour indices
        }
        statflag = 0L;

        birth = 0;
        if(gol[ij]) {                                                       // death/survival
            if (!((survivemask>>s)&1ull)) fprintf(stderr,"survival mask restricts birth for s= %d\n",s);
            if((((gene>>(s<<(ncoding-1))) & ncodingmask)==ncodingmask) && ((survivemask>>s)&1ull)) {       // survive
                statflag |= F_survival;
                newgol[ij]  = gol[ij];                                      // new game of life cell value same as old
                newgolg[ij] = golg[ij];                                     // gene stays same
            }
            else {                                                          // gene bit = 0 => death
                statflag |= F_death;
                newgol[ij]  = 0ull;                                         // new game of life cell value dead
                newgolg[ij] = 0ull;                                         // gene dies
                hashdeletegene(golg[ij],"step %d hash delete error 2 in update, gene %llx not stored\n");
            }
        }
        else{                                                               // birth/nobirth
            if(s==3)  {                                                     // birth
                birth = 1;                                                  // in this case, gene birth is coupled to gol state birth
                newgol[ij]  = 1ull;
                statflag |= F_birth;
            }
        }
        if(birth){                                                          // birth : we could in fact decouple the criteria for gene birth from that for gol state
            RAND128P(randnr);                                               // inline exp so compiler recognizes auto-vec,
            for(k=0;k<s;k++) livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];   // live neighbour genes
            // figure out newgene by random choice for now
            idx = ((randnr>>52)&0xf) % s;                                   // choose random for now
            newgene = livegenes[idx];
            // compute random events for single bit mutation, as well as mutation position nmut
            randnr2 = (randnr & pmutmask);                                  // extract bits from randnr for random trial for 0 on pmutmask
            r2 = randnr2?0ull:1ull;                                         // 1 if lowest nlog2pmut bits of randnr are zero, else zero
            nmut = (randnr >> 57) & 0x3f;                                   // choose mutation position for length 64 gene : from bits 57:62 of randnr
            // nmut = nmut & NGENE;          CHANGED
            // complete calculation of newgol and newgolg, including mutation
            ancestor = newgene;
            newgene = newgene ^ (r2<<nmut);                                 // introduce single mutation with probability pmut = probmut
            newgol[ij]  =  1L;                                              // new game of life cell value: alive
            newgolg[ij] =  newgene;                                         // if birth then newgene
            statflag = statflag | F_birth;
            if (r2) statflag = statflag | F_mutation;
            if(gol[ij]) hashdeletegene(golg[ij],"step %d hash delete error 1 in update_lut_sum, gene %llx not stored\n");
            hashaddgene(newgene,ancestor);
        }
        if(gol[ij]) statflag |= F_golstate;
        golgstats[ij] = statflag;
    }  // end for ij

    for (ij=0; ij<N2; ij++) {       // complete missing hash table records of extinction and activities
        if(gol[ij]) hashgeneextinction(golg[ij],"hash extinction storage error in update, gene %llx not stored\n");
        if(newgol[ij]) hashgeneactivity(newgolg[ij],"hash activity storage error in update, gene %llx not stored\n");
    }

    for (ij=0; ij<N2; ij++) {        // update lattices
        gol[ij] = newgol[ij];        // copy new gol config to old one
        golg[ij] = newgolg[ij];      // copy new genes to old genes
        golgstats[ij] = newgolgstats[ij]; // copy new golg stat flags to old genes
    }
}

//------------------------------------------------------- update -----------------------------------------------------------------------------------
void update(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]){
    /* update GoL for toroidal field which has side length which is a binary power of 2 */
    /* encode without if structures for optimal vector treatment */

    int s, s0, k, k1, nmut;
    unsigned int kch,checkmasks,checknongol,rulemod1,rulemod2,rulemod1ij,rulemod2ij;
    int nb[8], nbc, nbch, ij, i, j, jp1, jm1, ip1, im1;
    uint64_t g, gs, nb1i, nb2i, randnr, randnr2, r2;
    uint64_t nbmask, nbmaskr;
    uint64_t newgene, ancestor, livegenes[3];
    uint64_t s2or3, birth, statflag, nextgolstate;

    totsteps++;
    if(!(totsteps%10)) fprintf(stderr,"iteration step %d\r",totsteps);
    rulemod1= rulemod & 0x1ull;                                             // allow GoL rule modification
    rulemod2= (rulemod & 0x2)>>1;                                           // allow GoL rule modification masking
    checkmasks = rulemod2 && (repscheme & R_6_checkmasks);                  // enable genetic masking of neighbours
    checknongol = repscheme & R_7_nongolstat;                               // enable control of successive rule departures

    for (ij=0; ij<N2; ij++) {                                               // loop over all sites of 2D torus with side length N
	    i = ij & Nmask;  j = ij >> log2N;                                   // row & column
	    jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
	    ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                         // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;           // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (k=0,s=0,nb1i=0;k<8;k++) {                                      // packs non-zero nb indices in first up to 8*4 bits
            gs=gol[nb[k]];                                                  // whether neighbor is alive
            s += gs;                                                        // s is number of live nbs
            nb1i = (nb1i << (gs<<2)) + (gs*k);                              // nb1i is packed list of live neighbour indices
        }
        s0 = s;                                                             // record unmodified value of s for comparison
        s2or3 = (s>>2) ? 0ull : (s>>1);                                     // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
        nextgolstate = s2or3 ? (gol[ij] ? 1ull : (s&0x1ull ? 1ull : 0ull )) : 0ull;   // GoL standard calculation next state

        statflag = 0ull;
        rulemod2ij = rulemod2;
        rulemod1ij = rulemod1;
        nbmask = 0;
        if(s>1) {
          if(checknongol) {                                                 // check for non GoL changed states
            if(repscheme & R_8_nongolstatnbs) {                             // check in neighborhood
                int sng;                                                    // sum of nbs with non GoL rule change bit set
                for (k=0,sng=0;k<8;k++) {                                   // calc number of neighbours resulting from nongol rule change
                    sng += (golgstats[nb[k]]&F_nongolchg)?1:0;              // neighbors with state set by a non GoL change
                }
                if(sng) rulemod2ij = rulemod1ij = 0ull;
            }
            if(golgstats[ij]&F_nongolchg) rulemod2ij = rulemod1ij = 0ull;   // if central state the result of a non GoL rule change
          }
          if(checkmasks&&rulemod2ij) {                                      // recalculate effective new s as less than original sum sm
            for (k=0,nbmaskr=0;k<8;k++) {                                   // depending on rotated overlay of masks in live neighbour genes
                if(gol[nb[k]]) {
                    g= golg[nb[k]];                                         // fnal gene has ncoding 0s then 8 bit mask
                    //if(!((g>>8)&codingmask)) nbmaskr |= g&0xff;           // tried |=, &=, ^= .
                    if(!((g>>8)&codingmask)) nbmaskr |= (0x2<<(g&0x7L))&0xff;  // only one bit on for gene masks in this version: not self, so 0x2L

                }
                nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);           // 8 bit rotate right
            }
            //nbmaskr = 0ull;       // DEBUG check that default behaviour if genes no influence : PASSED
            for (k=0,s=0,nb1i=0ull;k<8;k++) {                               // recalculate sum and nb1i using combined mask
                gs =gol[nb[k]] & (((~nbmaskr)>>k)&0x1ull);                  // if mask bit set, count as if dead
                nbmask |= gs<<k;                                            // also calculate nbmask for use below
                s += gs;
                nb1i = (nb1i << (gs<<2)) + (gs*k);
            }
            if(s!=s0) s2or3 = (s>>2) ? 0ull : (s>>1);                       // redo s2or3 calculation for modified s
          }
        }
        if (s2or3) {                                                        // if 2 or 3 neighbours alive
            birth = 0ull;
            newgene = 0ull;
            if (s&0x1ull) {  // s==3                                        // allow birth (with possible overwrite)
              statflag |= F_3_live;                                         // record instance of 3 live nbs
              if ((0x1ull&overwritemask)|(0x1ull&~gol[ij]) ) {                  // central site empty or overwrite mode
                birth = 1ull;                                               // birth flag
                for(k=0;k<s;k++) livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]]; // live neighbour genes
                if (!(checkmasks&&rulemod2ij)) for (k=0,nbmask=0;k<8;k++) nbmask |= (gol[nb[k]]<<k);  // 8-bit mask of GoL states of 8 nbs, clockwise from top left
                statflag |= F_3_livenbs & (nbmask<<16);                     // record live neighbour pattern
                if((livegenes[0]^livegenes[1])|(livegenes[0]^livegenes[2])) { // genes not all same, need ancestor calculation
                  selectdifft3(nbmask, nb, gol, &kch);nbch=nb[kch&0x7];
                  if (repscheme & R_0_2sel_3live) {                         // execute selective replication of one of two otherwise unchosen live genes
                      nb2i = 0ull;
                      for(k1=k=0;k<3;k++) {                                 // choice of two other live genes for possible ancestor
                          nbc=(nb1i>>(k<<2))&0x7;
                          if(nb[nbc]!=nbch) nb2i = (nbc<<(k1++<<2))+nb2i;
                      }
                      if (repscheme & R_2_2ndnbs) selectone_nbs(s,nb2i,nb,gol,golg,&birth,&newgene);
                      else selectone(s,nb2i,nb,golg,&birth,&newgene);
                      if (birth==0ull) {                                    // optional reset of ancestor & birth if no ancestors chosen in selectone
                        if((repscheme & R_3_enforce3birth)||rulemod1ij)  {  // birth cannot fail or genes don't matter or no modification to gol rules
                            newgene = golg[nbch];
                            birth = 1ull;
                        }
                      }
                      else statflag |= F_2select;                           // ancestor has been chosen in selectone
                  }
                  else {
                      newgene = golg[nbch];
                  }
                  // if (newgene == 0ull) fprintf(stderr,"step %d Warning, new gene zero: nbmask %llx nbmaskrm %llx kmin %d gol %llx golg %llx newgene %llx ij %d\n",totsteps,nbmask,nbmaskrm,kmin,gol[nb[kmin]],golg[nb[kmin]],newgene,ij);
                } // end if not all live neighbors the same
                else {
                    statflag |= F_3g_same;
                    newgene = livegenes[0];                                 // genes all the same : copy first one
                    if((~repscheme&R_3_enforce3birth) && rulemod1ij) birth = 0ull;   // no birth for 3 identical genes
                }
              } // end central site empty or overwrite mode
            }  // end if s==3
            else {  // s==2                                                 // possible birth as exception to GoL rule
                statflag |= F_2_live;
                if (rulemod1ij||gol[ij]) {                                  // rule departure from GOL allowed or possible overwrite
                    if ((0x1ull&(overwritemask>>1))||!gol[ij]) {            // either overwrite on for s==2 or central site is empty
                        //for (k=0;k<s;k++) livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]]; // live gene at neighbour site
                        if (repscheme & R_2_2ndnbs) selectone_nbs(s,nb1i,nb,gol,golg,&birth,&newgene);
                        else if (repscheme & R_5_2uniquepos) {
                            nbmask = (0x1ull<<(nb1i&0x7)) + (0x1ull<<((nb1i>>4)&0x7));
                            selectdifft2(nbmask, nb, gol, &kch);
                            newgene = golg[nb[kch]];
                            if(repscheme & R_9_12_2birth_k1_4) {
                                if(repscheme & (R_9_2birth_k0<<(kch-1))) birth = 1ull;
                                else birth = 0ull;
                            }
                            else birth = 1ull;
                        }
                        else selectone(s,nb1i,nb,golg,&birth,&newgene);
                        if(!birth && (repscheme&R_4_enforce2birth)) {
                            nbmask = (0x1ull<<(nb1i&0x7)) + (0x1ull<<((nb1i>>4)&0x7));
                            selectdifft2(nbmask, nb, gol, &kch);
                            newgene = golg[nb[kch]];
                            if(repscheme & R_9_12_2birth_k1_4) {
                                if(repscheme & (R_9_2birth_k0<<(kch-1))) birth = 1ull;
                                else birth = 0ull;
                            }
                            else birth = 1ull;
                        }
                        if (birth) statflag |= F_2select;
                    }
                }
            }

            if(birth){
                // if (gol[ij]) fprintf(stderr,"birth overwrite event ij %d newgene %llu s %llu\n",ij,newgene,s);
                RAND128P(randnr);                                           // inline exp so compiler recognizes auto-vec,
                // compute random events for single bit mutation, as well as mutation position nmut
                randnr2 = (randnr & pmutmask);                // extract bits from randnr for random trial for 0 on pmutmask
                r2 = (!pmutmask||randnr2)?0ull:1ull;                        // 1 if lowest nlog2pmut bits of randnr are zero, else zero
                nmut = (randnr >> 56) & 0x3f;                               // choose mutation position for length 64 gene : from bits 56:61 of randnr
                // complete calculation of newgol and newgolg, including mutation
                ancestor = newgene;
                newgene = newgene ^ (r2<<nmut);                             // introduce single mutation with probability pmut = probmut
                // if(newgene^ancestor) fprintf(stderr,"newgene %llu written at step %d with ancestor %llu\n",newgene, totsteps,ancestor);     // DEBUG

                if(gol[ij]) {                                               // central old gene present: overwritten
                    hashdeletegene(golg[ij],"step %d hash delete error 1 in update, gene %llx not stored\n");
                }
                hashaddgene(newgene,ancestor);

                newgol[ij]  =  1ull;                                        // new game of life cell value: alive
                newgolg[ij] =  newgene;                                     // if birth then newgene
                statflag = statflag | F_birth;
                if (r2) statflag = statflag | F_mutation;
                // if(newgene==0ull) fprintf(stderr,"error in writing newgene, previous = %llx, statflag = %llx\n",golg[ij],statflag);
            } // end birth
            else {
                if ((survival&s&0x1ull)|((survival>>1)&(~s)&0x1ull)|((~rulemod1ij)&0x1ull)) {// (surv bit 0 and s==3) or (surv bit 1 and s==2) or not rulemod1ij
                // if ((survival&s&0x1ull)|((survival>>1)&(~s)&0x1ull)) {   // survival bit 0 and s==3, or (survival bit 1 and s==2)
                    newgol[ij]  = gol[ij];                                  // new game of life cell value same as old
                    newgolg[ij] = golg[ij];                                 // gene stays as before, live or not
                    if(gol[ij]) statflag |= F_survival;
                }
                else {
                    if(gol[ij]) {                                           // death : need to update hash table
                        hashdeletegene(golg[ij],"step %d hash delete error 2 in update, gene %llx not stored\n");
                    }
                    newgol[ij]  = 0ull;                                     // new game of life cell value dead
                    newgolg[ij] = 0ull;                                     // gene dies or stays dead
                    if(gol[ij]) statflag |= F_death;
                }
            } // end no birth
        }  // end if s2or3
        else {                                                              // else not birth or survival, 0 values for gol and gene
            if(gol[ij]) {                                                   // death : need to update hash table
                hashdeletegene(golg[ij],"step %d hash delete error 3 in update, gene %llx not stored\n");
            }
	        newgol[ij]  = 0ull;                                             // new game of life cell value
	        newgolg[ij] = 0ull;                                             // gene dies
            if(gol[ij]) statflag |= F_death;
        }
        if(gol[ij]) statflag |= F_golstate;
        if(newgol[ij]^nextgolstate) statflag |= F_notgolrul;
        if(gol[ij]^newgol[ij]) {
            statflag |= F_golchange;
            if(statflag&F_notgolrul) statflag |= F_nongolchg;
        }
        else if (golgstats[ij]&F_nongolchg) statflag |= F_nongolchg;        // maintain non-GoL chg status until state changed by GoL rule
        newgolgstats[ij] = statflag;
    }  // end for ij

    for (ij=0; ij<N2; ij++) {       // complete missing hash table records of extinction and activities
        if(gol[ij]) hashgeneextinction(golg[ij],"hash extinction storage error in update, gene %llx not stored\n");
        if(newgol[ij]) hashgeneactivity(newgolg[ij],"hash activity storage error in update, gene %llx not stored\n");
    }

    for (ij=0; ij<N2; ij++) {        // update lattices
	    gol[ij] = newgol[ij];        // copy new gol config to old one
	    golg[ij] = newgolg[ij];      // copy new genes to old genes
        golgstats[ij] = newgolgstats[ij]; // copy new golg stat flags to old genes
    }

}

//------------------------------------------------------- stats routines -----------------------------------------------------------------------------------
void tracestats(uint64_t gol[],uint64_t golg[], uint64_t golgstats[], int NN2) { // trace various stats over time of the simulation
    int ij,cnt,k,d,dc,gt[4],st[10];
    uint64_t gene,statflag;
  
    if (statcnts == arraysize) {                                            // relallocate memory for arrays : double size
        arraysize*=2;
        livesites = (int *)realloc(livesites, arraysize * sizeof(int));
        genestats = (int *)realloc(genestats, arraysize * 4 * sizeof(int));
        stepstats = (int *)realloc(stepstats, arraysize * 10 * sizeof(int));
        if (nhistG==nstatG) configstats = (int *)realloc(configstats, arraysize * Noff * sizeof(int));
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
            switch (selection) {
                case 0: if(d==64) d--; dc=(d>>4)&0x3;break;
                case 1: if(d==64) d--; dc=(d>>4)&0x3;break;
                case 2:
                case 3: dc=d&0x3;break;
                default:if(d==64) d--; dc=(d>>4)&0x3;
            }
            gt[dc]++;
        }
    }

    livesites[statcnts] = cnt;
    for(d=0;d<4;d++) genestats[statcnts*4+d]=gt[d];
    for(k=0;k<10;k++) stepstats[statcnts*10+k]=st[k];
    if (nhistG==nstatG) for(k=0;k<10;k++) configstats[statcnts*Noff+k]=histo[k];
    statcnts++;
}

void get_stats(int outstats[], int outgtypes[], int outstepstats[], int outconfigstats[], int numStats ){
    int i;
    if(numStats > arraysize){
        fprintf(stderr,"Ack! numStats = %d  > arraysize = %d\n",numStats,arraysize);
        exit(1);
    }
    for(i=0; i<numStats; i++) outstats[i] = livesites[i];
    for(i=0; i<4*numStats; i++) outgtypes[i] = genestats[i];
    for(i=0; i<10*numStats; i++) outstepstats[i] = stepstats[i];
    if (nhistG==nstatG) for(i=0; i<Noff*numStats; i++) outconfigstats[i] = configstats[i];
}

//------------------------------------------------------- countconfigs -----------------------------------------------------------------------------------
void countconfigs(){        // count translated configs specified by offset array
    // for non zero nb patterns we count if the same pattern occurs at the central site and its offset address in space time
    int ij,i,j,k,t,ip1,im1,jp1,jm1,o,ij1,i1,j1;
    uint64_t nbmask1,nbmask2;
    int nb[9];
    uint64_t *pl, *pl0;

    for (o=0;o<Noff;o++) histo[o] = 0;

    pl0 = planes[curPlane];
    for (ij=0; ij<N2; ij++) {                                               // loop over all sites of 2D torus with side length N
        i = ij & Nmask;  j = ij >> log2N;                                   // row & column
        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                         // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;           // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        nb[8]=ij;
        for (k=0,nbmask1=0ull;k<9;k++) nbmask1 |= (*(pl0+nb[k]))<<k;            // calculate nbmask for comparison below
        if(nbmask1) {
            for (o=0;o<Noff;o++) {
                i1 = (i+offsets[o][0]) & Nmask;                              // periodic boundary conditions for N power of 2
                j1 = (j+offsets[o][1]) & Nmask;
                t = (curPlane+offsets[o][2]) & (numPlane-1);                // periodic in numPlane if this is power of 2
                pl = planes[t];
                ij1 = i1 + j1*N;
                jp1 = ((j1+1) & Nmask)*N; jm1 = ((j1-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
                ip1 =  (i1+1) & Nmask; im1 =  (i1-1) & Nmask;                         // toroidal i+1, i-1
                nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;           // new order of nbs
                nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
                nb[8]=ij1;
                for (k=0,nbmask2=0ull;k<9;k++) nbmask2 |= (*(pl+nb[k]))<<k;                          // also calculate nbmask for use below
                if(nbmask1==nbmask2) histo[o]++;
            }
        }
    }
}

//------------------------------------------------------- get_ ... -----------------------------------------------------------------------------------
void get_histo(int outhisto[],int numHistoC){
    int i;
    if(numHistoC != numHisto){
        fprintf(stderr,"Ack! numHisto = %d  != numHistoC = %d\n",numHisto,numHistoC);
        exit(1);
    }
    for(i=0; i<numHisto; i++) outhisto[i] = histo[i];
}

#ifdef HASH
void get_activities(uint64_t actgenes[],int activities[],int ngenesp[]) {
    int k, nlivegenes, nspecies;

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );
    fprintf(stderr,"The number of different species that have ever existed is %d\n",nspecies);
    
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
#else
void get_activities(uint64_t actgenes[],int activities[],int ngenesp[]) {}
#endif

//------------------------------------------------------- genelife_update -----------------------------------------------------------------------------------
void genelife_update (int nsteps, int nhist, int nstat) {
    /* update GoL for toroidal field which has side length which is a binary power of 2 */
    /* encode without if structures for optimal vector treatment */
    int t;
    uint64_t *newgol, *newgolg;
    int *activities,*gindices,*popln,*birthsteps,nspecies,ngenealogydeep;
    uint64_t *genes;
    int activitieshash(int gindices[], uint64_t genes[], int popln[], int activities[],int col);   /* count activities of all currently active species */
    int genealogies(int gindices[], uint64_t genes[], int popln[], int activities[], int birthsteps[]);   /* genealogies of all currently active species */
    
    nhistG = nhist;
    nstatG = nstat;
    for (t=0; t<nsteps; t++) {
        newgol = planes[newPlane];
        newgolg = planesg[newPlane];

        if (selection<10)       update(gol,golg,newgol,newgolg);              // calculate next iteration with selection
        else if (selection <12) update_gol16(gol,golg,newgol,newgolg);        // calculate next iteration for 2-16x multiplane version
        else if (selection==12) update_gol64(gol,golg,newgol,newgolg);        // calculate next iteration for 64x multiplane version
        else if (selection==13) update_gol2(gol,golg,newgol,newgolg);         // calculate next iteration for 2xgol multiplane version
        else if (selection==14) update_lut_sum(gol,golg,newgol,newgolg);      // calculate next iteration for lut sum version
        if(nhist && (totsteps%nhist == 0)) countconfigs();                    // count configurations
        if(nstat && (totsteps%nstat == 0)) tracestats(gol,golg,golgstats,N2); // time trace point

        gindices=NULL;activities=NULL;genes=NULL;popln=NULL;birthsteps=NULL;  // these arrays are mallocated in activitieshash or genealogies
        nspecies=activitieshash(gindices, genes, popln, activities,1);        // colors acttrace and returns current population arrays
        // possible further use of returned current gene population data here
        if(nspecies<0) fprintf(stderr,"error returned from activitieshash\n");
        free(gindices);free(activities);free(genes);free(popln);free(birthsteps);// free arrays after use

        gindices=NULL;activities=NULL;genes=NULL;popln=NULL;birthsteps=NULL;  // these arrays are mallocated in activitieshash or genealogies
        ngenealogydeep=genealogies(gindices, genes, popln, activities,birthsteps); // colors genealogytrace
        // possible further use of returned current gene population data here
        if(ngenealogydeep<0) fprintf(stderr,"error returned from genealogies\n");
        free(gindices);free(activities);free(genes);free(popln);free(birthsteps);// free arrays after use

        curPlane = (curPlane +1) % numPlane;            // update plane pointers to next cyclic position
        newPlane = (newPlane +1) % numPlane;
        gol = planes[curPlane];                         // get planes of gol,golg data
        golg = planesg[curPlane];
        totdisp++;                                      // currently every step is counted for display in activities
    }
}

//------------------------------------------------------- initialize_planes ---------------------------------------------------------------------------
void initialize_planes(int offs[],  int No) {
    int i,j,idx;
    static int notfirst = 0;

    curPlane = 0;
    newPlane = 1;
    if (notfirst)   return;     // need to fix memory free at two levels unless this fix
    notfirst = 1;
    
    if(No%3 !=0) fprintf(stderr,"Ouch!  Size of offsets array not a multiple of 3 as expected.");
    Noff = No/3;		// Noff global
    if(Noff>24){
        fprintf(stderr,"Too many offsets!  Max=24");
        exit(1);
    }
    numHisto = Noff;
    histo = (int *) calloc(numHisto,sizeof(int));

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
    if(tmx>0)   {                // exit if positive values of z in (x,y,z) offsets have been entered : these look into future
        fprintf(stderr,"Error: offsets looking into future not allowed ------- tmx = %d, tmn = %d",tmx,tmn);
        exit(1);
    }    tall = tmx-tmn;
    // numPlane = 2 + tall;    // numPlane >= 2
    if(tall>4) numPlane = 8;                            // only numPlane values which are powers of two are allowed
    else if (tall > 2) numPlane = 4;
    else numPlane = 2;
    if (numPlane>maxPlane) {
        fprintf(stderr,"Not enough planes defined by maxPlane for given offsets (need > %d)\n",numPlane);
        exit(1);
    }

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
}

//------------------------------------------------------- readFile and writeFile ---------------------------------------------------------------
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

//------------------------------------------------------- initialize ---------------------------------------------------------------------------
void initialize(int runparams[], int nrunparams, int simparams[], int nsimparams) {
#ifdef HASH
    int hcnt;
#endif
    int ij,ij1,i0,j0,i,j,Nf,k,cnt,icf,nstartgenes;
    unsigned int np;
    uint64_t g,pmask;

    static unsigned int rmask = (1 << 15) - 1;
    static int notfirst = 0;
    // Range: rand returns numbers in the range of [0, RAND_MAX ), and RAND_MAX is specified with a minimum value of 32,767. i.e. 15 bit

    uint64_t startgenes[16];
    char *golgin;
    
    srand(1234567);
    state[0] = rand();state[1] = rand();
    cnt = 0;
    totsteps = 0;
    totdisp = 0;
    statcnts = 0;
    
    // writeFile("genepat.dat");

    rulemod = runparams[0];
    repscheme = runparams[1];
    selection = runparams[2];
    overwritemask = runparams[3];
    survival = runparams[4];
    colorfunction = runparams[5];
    inittype = runparams[6];

    nlog2pmut = simparams[0];
    if(nlog2pmut>56) nlog2pmut=0;                           // need to use top 6-8 bits of 64 bit random nr for position
    pmutmask = (0x1ull << nlog2pmut) - 0x1ull;              // NB if nlogpmut==0, pmutmask==zero, no mutation.
    initial1density = simparams[1];
    initialrdensity = simparams[2];
    ncoding = simparams[3];
    codingmask = (0x1ull<<ncoding)-0x1ull;
    startgenechoice = simparams[4];
    
    fprintf(stderr,"runparams %d %d %d %d %d %d %d\n",runparams[0],runparams[1],runparams[2],
                                         runparams[3],runparams[4],runparams[5],runparams[6]);
    fprintf(stderr,"simparams %d %d %d %d %d\n",simparams[0],simparams[1],simparams[2],simparams[3],simparams[4]);
    fprintf(stderr,"pmutmask %llx (NB 0 means no mutation)\n",pmutmask);
                                                                              // not yet useable with selection==12
    np = (repscheme && repscheme<=16) ? repscheme : 16;                       // number of planes used, 16 if 0
    pmask = (np==16) ? 0xffffffffffffffffull : (1ull<<(np<<2))-1;                // mask for bits included in np planes
    if(selection<10 || selection > 11) np = 64;
    if (selection == 10) gene0=0xfedcba9876543210ull&pmask;
    else                 gene0=0ull;

    nstartgenes = 8;
    if (selection >= 10 && selection <=12) nstartgenes = 16;
    switch (selection) {
        case 0: for (k=0;k<4;k++) {startgenes[k]=0xf0f0f0f0f0f0f0f0;startgenes[k+4]=0x0f0f0f0f0f0f0f0f;} break;
        case 1: for (k=0;k<8;k++) startgenes[k]=((0x1ull<<k*3)-1ull)<<20;break;
        case 2:
        case 3: for (k=0;k<8;k++) startgenes[k]=(((0x1ull<<20)-1ull)<<20)+((0x1ull<<k)-0x1ull);break;
        case 4: for (k=0;k<8;k++) {g = 0xfffff0ull + k; startgenes[k] = k<4 ? g : ~g;} break;
        case 5: for (k=0;k<8;k++) {g = 0xf0ull + k; startgenes[k]= k<4 ? g : (~g)|(0xfull<<16);} break;
        case 10:for (k=0;k<16;k++) startgenes[k] = gene0;                     // first set all startgenes as uncoupled
                for (k=0;k<np;k++) startgenes[k] = (startgenes[k] & (~(0xfull<<(k<<2)))) | ((k+1)%np)<<(k<<2);break; // couple plane k+1 to k in startgene k
        case 11:for (k=0;k<16;k++) startgenes[k] = gene0;                     // first set all startgenes as uncoupled
                for (k=0;k<np;k++) startgenes[k] = 6<<(k<<2); break;          // coupled to neighbour plane before and after
        case 12:for (k=0;k<16;k++) startgenes[k] = 0x3ull<<(k<<2);break;
        case 13:for (k=0;k<8;k++) startgenes[k]=((0x1ull<<k*3)-1ull)<<20;break;
        case 14:for (k=0;k<8;k++) startgenes[k]=0xff000000+0x0c;
        case 6:
        case 7:
        default: for (k=0;k<8;k++) startgenes[k]=(0x1ull<<(4+k*8))-1ull;
    }

    if ( livesites !=NULL) {
        free(livesites);
        livesites = NULL;
        free(genestats);
        genestats = NULL;
        free(stepstats);
        stepstats = NULL;
        if ( configstats != NULL) {
            free(configstats);
            configstats = NULL;
        }
    }
    arraysize = startarraysize;
    livesites = (int *) malloc(arraysize * sizeof(int));
    genestats = (int *) malloc(arraysize * 4 * sizeof(int));
    stepstats = (int *) malloc(arraysize * 10 * sizeof(int));
    if (nhistG==nstatG) configstats = (int *) malloc(arraysize * Noff * sizeof(int));
    curPlane = 0;                                           // if we rerun initialize, we want to restart plane cycling from zero
    newPlane = 1;
    gol = planes[curPlane];
    golg = planesg[curPlane];

#ifdef HASH
    if(notfirst) hashtable_term(&genetable);
    hashtable_init(&genetable,sizeof(genedata),N2<<2,0);     // initialize dictionary for genes
#endif
    notfirst = 1;
    if (inittype==1) {           // input from file genepat.dat with max size of 32*32 characters
        golgin = (char *) malloc(32* 32 * sizeof(char));
        icf=readFile(golgin, "genepat.dat");
        if (icf != 32*32) {
            icf = 0;
            fprintf(stderr,"error reading file, %d not 32*32 chars\n",icf);
        }
        for (ij=0; ij<N2; ij++) {
            gol[ij] = 0;
            golg[ij] = gene0;
        }
        for (ij1=0; ij1<32*32; ij1++) {
            if(N<32) {fprintf(stderr,"Error, array dim %d too small for file array dim %d\n",N,32);break;}
            ij=(N>>1)-16+(ij1&0x1f)+ N*((N>>1)-16+(ij1>>5));
            if (golgin[ij1] > 0)    {                   // if live cell
                gol[ij] = 1ull;
                if(golgin[ij1] <= 8 ) golg[ij] = startgenes[golgin[ij1]-1];
                else if (golgin[ij1]>='0' && golgin[ij1]<'8') golg[ij] = startgenes[golgin[ij1]-'0'];
                else golg[ij] = startgenes[7];
                cnt++;
                golgstats[ij] = 0;
            }
            else {
                gol[ij] = gene0;
                golgstats[ij] = 0;
            }
            // if (golg[ij] == 0 && gol[ij] != 0) fprintf(stderr,"zero gene at %d\n",ij);
        }

    }
    else {  // inittype gives linear size of random block for initialization (0 => full frame, as before)
        Nf = inittype;
        if (Nf==0 || Nf>N) Nf=N;
        for (ij=0; ij<N2; ij++) {
            gol[ij] = 0ull;
            golg[ij] = gene0;
            golgstats[ij] = 0;
        }
        i0 = j0 = (N>>1)-(Nf>>1);
        for (i=0; i<Nf; i++) {
            for (j=0; j<Nf; j++) {
                ij=i0+i+N*(j0+j);
                if((selection<10)||(selection==14)) gol[ij] = ((rand() & rmask) < initial1density)?1ull:0ull;
                else if(selection<13) for (k=0;k<np;k++) gol[ij] |= ((rand() & rmask) < initial1density)?(1ull<<(k<<2)):0ull;
                else if (selection==13) for (k=0;k<2;k++) gol[ij] |= ((rand() & rmask) < initial1density)?(1ull<<k):0ull;
            }
        }
        for (ij=0; ij<N2; ij++) {
            g = 0ull;
            if (gol[ij]||(selection>=10)) { // if live cell or multiplane, fill with random genome g or randomly chosen startgene depending on initialrdensity
                if (((unsigned) rand() & rmask) < initialrdensity) {for (k=0; k<((np>16)?np:(np<<2)); k++) g = (g << 1) | (rand() & 0x1);g=gene0^g;}  // all 64 for 1selection<10 or selection>11 or np=16
                else if (startgenechoice == nstartgenes) g = startgenes[0xf & rand() & (nstartgenes-1)];
                else if (startgenechoice > nstartgenes) fprintf(stderr,"startgenechoice %d out of range\n",startgenechoice);
                else g = startgenes[0xf & startgenechoice & (nstartgenes-1)];
                cnt++;
            }
            golg[ij] = g;
        }
    }
    for (ij=0; ij<N2; ij++) {
        acttrace[ij]=rootgene;                 // initialize activity traces to root gene
        genealogytrace[ij] = rootgene;
    }

    for (ij=0; ij<N2; ij++) {
        if(gol[ij]||(selection==10)||(selection==11)||((gol[ij]>>1)&&(selection==13))) {
            hashaddgene(golg[ij],rootgene);
        }
    }

#ifdef HASH
    // it is possible to enumerate keys and values
    hcnt=hashtable_count(&genetable);
    genotypes = hashtable_keys( &genetable );
    fprintf(stderr,"population size %d with %d different genes\n",cnt,hcnt);
#endif
}

//------------------------------------------------------- set ...---------------------------------------------------------------------------
void set_colorfunction(int colorfunctionval) {
    if(colorfunction>7) fprintf(stderr,"error colorfunction value passed %d too large\n",colorfunctionval);
    else     colorfunction = colorfunctionval;
}

int setget_act_ymax(int actymax) {                  // sets ymax for activities only if argument nonzero, reads old value
    int ymaxold;
    ymaxold = ymax;
    ymax = actymax;
    return(ymaxold);
}

//------------------------------------------------------- get ... ---------------------------------------------------------------------------
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

void get_acttrace(uint64_t outgolg[], int NN){
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolg[ij] = acttrace[ij];
    }
}

void get_genealogytrace(uint64_t outgolg[], int NN){
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolg[ij] = genealogytrace[ij];
    }
}

int get_nspecies(){
    int k,nspecies,nspeciesnow;
    nspecies = hashtable_count(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );

    for (k=0,nspeciesnow=0; k<nspecies; k++)
        if(geneitems[k].popcount) nspeciesnow++;
    return(nspeciesnow);
}


void get_curgolgstats(uint64_t outgolgstats[], int NN){
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolgstats[ij] = golgstats[ij];                       // Note that golgstats is not dealt with in planes !
    }
}

//------------------------------------------------------- countspecies ---------------------------------------------------------------------------
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

void countspecies1(uint64_t gol[], uint64_t golg[], int N2) {     /* counts numbers of all different species using qsort first */
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
    fprintf(stderr,"The number of different species (countspecies) is %d\n",nspecies);

    for (k=0,ij=0;k<nspecies;k++) {     // now condense array to give only different genes with counts
        // printf("species %4d with gene %x has counts %d\n",k, golgs[ij],counts[k]);
        golgs[k]=golgs[ij];
        ij = ij + counts[k];
    }

    for (k=0; k<nspecies; k++) { golgsc[k][0] = golgs[k];  golgsc[k][1] = counts[k];}  // initialize joint gene & count array
    qsort(golgsc, nspecies, sizeof(golgsc[0]), cmpfunc1);                   // sort in decreasing count order

#ifdef HASH
    for (k=1; k<nspecies; k++) {                            // check consistency of hash table data, assuming empty site gene is most frequent
        if((genedataptr = (genedata *) hashtable_find(&genetable, golgsc[k][0])) != NULL) {
                    if(genedataptr->popcount != golgsc[k][1])
                        fprintf(stderr,"popcount %llu <> %d hash error at k = %d\n",golgsc[k][1],genedataptr->popcount,k);
        }
        else fprintf(stderr,"countspecies popcount error, no entry in hash table\n");
    }
#endif

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
        else if ((selection == 2)||(selection == 3)) {                      // cyclic 4 species model
	        fitness = nones&0x3;                                            // fitness is species class
        }

        fprintf(stderr,"count species %d with gene %llx has counts %llu and %d ones, fitness %llu\n",k, golgsc[k][0],golgsc[k][1],nones,fitness);
    }
    fprintf(stderr,"rulemod\trepscheme\tselection\toverwritemask\tsurvival\n");
    fprintf(stderr,"%d\t%d\t\t%d\t\t%d\t\t%d\n",rulemod,repscheme,selection,overwritemask,survival);
    fprintf(stderr,"nlog2pmut\tinit1\tinitr\tncoding\tstartchoice\n");
    fprintf(stderr,"%d\t\t%d\t%d\t%d\t%d\n",nlog2pmut,initial1density,initialrdensity,ncoding,startgenechoice);
}

void countspecies() {                                                       // counts current species without using hash tables
    countspecies1(gol, golg, N2);
}

#ifdef HASH
int cmpfunc2 (const void * pa, const void * pb)
{
   return ( genotypes[*(int*)pa] > genotypes[*(int*)pb] ? 1 : -1);
}

int cmpfunc3 (const void * pa, const void * pb)
{
   return ( geneitems[*(int*)pa].popcount < geneitems[*(int*)pb].popcount ? 1 : -1);
}

void countspecieshash() {  /* counts numbers of all different species using qsort first */
    int k, *golgs, nspecies, nspeciesnow, nones;
    uint64_t last, fitness;
//    static int debugcnt = 0;

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );
    
//    fprintf(stderr,"Debug cnt %d in countspecieshash\n",debugcnt++);
    golgs = (int *) malloc(nspecies*sizeof(int));
    for (k=0; k<nspecies; k++) golgs[k] = k;  // initialize sorted genotype array to same order as hash table

    // qsort(golgs, nspecies, sizeof(int), cmpfunc2);                       // sort in increasing gene order
    qsort(golgs, nspecies, sizeof(int), cmpfunc3);                          // sort in decreasing count order

    // for (k=0; k<nspecies; k++) fprintf(stderr,"in countspecieshash genotype %d is %llx\n", k, genotypes[k]);
    for (k=0,nspeciesnow=0; k<nspecies; k++) {
        last = genotypes[golgs[k]];
        POPCOUNT64C(last, nones);
        fitness = 999L;
        if (selection == 0) {                                               // 2-live neighbor fitness is integer value
            fitness = last;
        }
        else if (selection == 1) {                                          // 2-live neighbor fitness is number of ones
            fitness = (uint64_t) nones;
        }
        else if ((selection == 2)||(selection == 3)) {                      // cyclic 4 species model
            fitness = nones&0x3;                                            // fitness is species class
        }

        if((genedataptr = (genedata *) hashtable_find(&genetable, last)) != NULL) {
            if(genedataptr->popcount) {
                fprintf(stderr,"count species %7d with gene %16llx has counts %7d and %3d ones, fitness %llu\n",k,last,
                    genedataptr->popcount,nones,fitness);
                nspeciesnow++;
            }
        }
        else {
            fprintf(stderr,"countspecieshash popcount error, no entry in hash table\n");
            fprintf(stderr,"count species %d with gene %llx has counts ?? and %d ones, fitness %llu\n",k,last,
                    nones,fitness);
        }
    }
    fprintf(stderr,"Iteration %d : %d different species of %d ever existed.\n",totsteps,nspeciesnow,nspecies);
    fprintf(stderr,"_________________________________________________________________\n");
    free(golgs);
    golgs = NULL;
}
#else
void countspecieshash() {}  /* dummy routine, no effect */
#endif

//------------------------------------------------------- activitieshash ---------------------------------------------------------------------------
#ifdef HASH
int activitieshash(int gindices[], uint64_t genes[], int popln[], int activities[], int col) {  /* count activities of all currently active species */
    int i, j, ij, ij1, x, nspecies, nspeciesnow;
    const int maxact = 10000;
    // int ymax1;
    uint64_t gene;
    
    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );

    if(gindices != NULL && col) free(gindices);
    for (i=0,nspeciesnow=0; i<nspecies; i++)
        if(geneitems[i].popcount) nspeciesnow++;
    if (col) gindices = (int *) malloc(nspeciesnow*sizeof(int));
    else if (nspeciesnow>10000) return(-1);                        // exit with error need to allocate more space in python
    for (i=j=0; i<nspecies; i++) {
        if(geneitems[i].popcount) {
            gindices[j]=i;                                         // if col is 0 then the array gindices must be passed with sufficient length
            j++;
        }
    }
    if (nspeciesnow > maxact || !col) {                             //sort with in order of decreasing population
        qsort(gindices, nspeciesnow, sizeof(int), cmpfunc3);        // sort in decreasing count order
    }
    if (nspeciesnow > maxact) nspeciesnow = maxact;

    if (col) {
        genes = (uint64_t *) malloc(nspeciesnow*sizeof(uint64_t));
        popln = (int *) malloc(nspeciesnow*sizeof(int));
        activities = (int *) malloc(nspeciesnow*sizeof(int));
    }

    for (i=0; i<nspeciesnow; i++) {
        genes[i]=genotypes[gindices[i]];
        popln[i]=geneitems[gindices[i]].popcount;
        activities[i]=geneitems[gindices[i]].activity;
    }
    
    if(!col) return(nspeciesnow);                                   // exit here unless doing display
    
    if (totdisp>=N) {                                               // 1 pixel to left scroll when full
        for(ij=0;ij<N2;ij++) {
            ij1 = ((ij+1)&Nmask)+((ij>>log2N)<<log2N);              // (i+1)%N+j*N;
            if(ij1>=N2) fprintf(stderr,"error in scroll of acttrace\n");
            acttrace[ij]=acttrace[ij1];
        }
        x=N-1;
    }
    else x=totdisp;

    for(i=0;i<N;i++) acttrace[x+i*N]=rootgene;                  // set column gray
    //for(i=ymax1=0;i<nspeciesnow;i++)
    //    ymax1 = activities[i]>ymax1 ? activities[i] : ymax1;
    // if (ymax1>ymax) ymax = ymax*2;     // autoscale of activities
    // if (ymax1<ymax/2) ymax = ymax/2;   // autoscale of activities
    for(i=0;i<nspeciesnow;i++) {
        activities[i] = N - (activities[i] * N) / ymax;
        gene = genes[i];
        ij = (x&Nmask)+activities[i]*N;
        if(ij > 0 && ij<N2)
            acttrace[ij] = gene;
        //else
        //     fprintf(stderr,"activity out of range\n");
    }
    return(nspeciesnow);
}
#else
int activitieshash(int gindices[], uint64_t genes[], int popln[], int activities[]) {return(0)}  /* dummy routine, no effect */
#endif

int get_sorted_popln_act( int gindices[], uint64_t genes[], int popln[], int activities[]) {
    int nspecies;
    nspecies=activitieshash(gindices, genes, popln, activities, 0);        // sets acttrace and returns current population arrays
    return(nspecies);
}

//------------------------------------------------------- genealogies ---------------------------------------------------------------------------
#ifdef HASH
int cmpfunc4 (const void * pa, const void * pb)
{
   return ( geneitems[*(int*)pa].firstbirthframe > geneitems[*(int*)pb].firstbirthframe ? 1 : -1);
}

int cmpfunc5 (const void * pa, const void * pb)                 // sort according to ancestry in genealogytrace
{
   int i1,i2,ij1,ij2,j;
   uint64_t gene1,gene2;
   i1=*(int*)pa; i2=*(int*)pb;
   for (j=0;j<genealogydepth;j++) {
        ij1 = i1+j*N; ij2 = i2+j*N;
        gene1=genealogytrace[ij1]; gene2=genealogytrace[ij2];
        if(gene1!=gene2) return((((gene1 > gene2) && (gene1!=rootgene)) || (gene2==rootgene)) ? 1 : -1);
    }
    return(0);
}


int genealogies(int gindices[], uint64_t genes[], int popln[], int activities[], int birthsteps[]) {  /* genealogies of all currently active species */
    int j, jmax, i, ij, nspecies, nspeciesnow, birthstep, *gorder;
    int j1, j2, j3, activity;
    uint64_t gene, ancgene, nextgene, *genealogytrace1;
    
    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );

    if(gindices != NULL) free(gindices);
    gindices = (int *) malloc(nspecies*sizeof(int));
    for (i=nspeciesnow=0; i<nspecies; i++)
        if(geneitems[i].popcount) nspeciesnow++;

    for (i=j=0; i<nspecies; i++) {
        if(geneitems[i].popcount) {
            gindices[j]=i;
            j++;
        }
        else gindices[nspeciesnow+i-j]=i;
    }
    // qsort(gindices, nspeciesnow, sizeof(int), cmpfunc4);// sort in increasing birthstep order
    qsort(gindices, nspeciesnow, sizeof(int), cmpfunc3);// sort in decreasing population size order
    
    if (nspeciesnow>N) nspeciesnow=N;               // can only display at most N species, chose oldest

    genes = (uint64_t *) malloc(nspeciesnow*sizeof(uint64_t));
    popln = (int *) malloc(nspeciesnow*sizeof(int));
    activities = (int *) malloc(nspeciesnow*sizeof(int));
    activitymax=0;
    birthsteps = (int *) malloc(nspeciesnow*sizeof(int));

    for (i=0; i<nspeciesnow; i++) {
        genes[i]=genotypes[gindices[i]];
        popln[i]=geneitems[gindices[i]].popcount;
        activities[i]=activity=geneitems[gindices[i]].activity;
        birthsteps[i]=geneitems[gindices[i]].firstbirthframe;
    }
    
    for(ij=0;ij<N2;ij++) genealogytrace[ij]=rootgene;             // set field to rootgene black
    ancgene=rootgene;                                             // never really used, but included to avoid unitialized warning
    birthstep=0;
    for (i=jmax=0; i<nspeciesnow; i++) {
        //j1=0;
        for (j=0;j<N;j++) {  // go back at most N links in genealogy
            if(j) {
                gene=ancgene;
                if(gene==rootgene) {
                    ancgene = rootgene;
                    birthstep = 0;
                    activity = 0;
                }
                else {
                    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
                        ancgene=genedataptr->firstancestor;
                        birthstep = genedataptr->firstbirthframe;
                        activity = genedataptr->activity;
                        if(activity>activitymax) activitymax=activity;
                    }
                    else fprintf(stderr,"ancestor not found in genealogies\n");
                }
            }
            else  {
                gene=genes[i];
                ancgene=geneitems[gindices[i]].firstancestor;
                birthstep=geneitems[gindices[i]].firstbirthframe;
                activity=geneitems[gindices[i]].activity;
                if(activity>activitymax) activitymax=activity;
            }
            if (gene == rootgene) break;                            // reached root, exit j loop
            ij = i+j*N;
            genealogytrace[ij]=gene;
        }
        if (j>jmax) jmax=j;
    }
    genealogydepth = jmax;
    
                                                                    //reverse ancestries to allow comparison at same number of speciations
    for (i=0; i<nspeciesnow; i++) {
        for(j=0;j<N;j++) {
            if (genealogytrace[i+j*N]==rootgene) break;
        }
        for(j1=0;j1<(j>>1);j1++) {
            gene=genealogytrace[i+(j-j1-1)*N];
            genealogytrace[i+(j-j1-1)*N]=genealogytrace[i+j1*N];
            genealogytrace[i+j1*N]=gene;
        }
    }
    gorder = (int *) malloc(N*sizeof(int));
    for (i=0; i<N; i++) gorder[i]=i;
    qsort(gorder, nspeciesnow, sizeof(int), cmpfunc5);              // sort according to ancestral lines
    for (i=0;i<N;i++) if((gorder[i]<0)||(gorder[i]>=N)) fprintf(stderr,"step %d error in gorder out of bounds at i = %d with value %d\n",totsteps,i,gorder[i]);

    genealogytrace1 = (uint64_t *) malloc(N*jmax*sizeof(uint64_t)); // allocate copy array of genetrace to allow sorting and modifications
    for(ij=0;ij<N*jmax;ij++) genealogytrace1[ij]=genealogytrace[ij];// copy active portion of genealogytrace to new array genealogytrace1
    for(ij=0;ij<N2;ij++) genealogytrace[ij]=rootgene;               // reinitialize genealogytrace to root gene before redrawing part of it

    if(colorfunction>=6) {                                          // time trace of genealogies
      for(i=0;i<nspeciesnow;i++) {
        for(j=0,j1=0;j<jmax;j++) {
            if(gorder[i]>=nspeciesnow) fprintf(stderr,"error in genealogies gorder at i=%d, order value %d out of range\n",i,gorder[i]);
            ij = gorder[i]+j*N;
            gene = genealogytrace1[ij];
            ij+=N;
            if(ij<N2) {
                nextgene = genealogytrace1[ij];
                if(nextgene==rootgene) birthstep=totsteps;
                else {
                    if((genedataptr = (genedata *) hashtable_find(&genetable, nextgene)) != NULL) birthstep = genedataptr->firstbirthframe;
                    else fprintf(stderr,"ancestor %llx not found at (%d,%d) in genealogies during birthstep extraction\n",nextgene,ij&Nmask,ij>>log2N);
                }
            }
            else birthstep=totsteps;
            j2 = birthstep*N/totsteps;
            for (j3=j1;j3<j2;j3++) {
                ij = i+j3*N;
                genealogytrace[ij]=gene;
            }
            j1 = j2;
        }
      }
      for(i=nspeciesnow;i<N;i++) for(j=0;j<N;j++) genealogytrace[gorder[i]+j*N]=rootgene;
    }
    else {                                                          // species changes only trace (colorfunction == 5)
      for(i=0;i<nspeciesnow;i++) {
        for(j=0;j<jmax;j++) {
            ij=i+j*N;
            genealogytrace[ij]=genealogytrace1[gorder[i]+j*N];
        }
      }
    }
    free(gorder);free(genealogytrace1);

    return(jmax);
}
#else
int genealogies(int gindices[], uint64_t genes[], int popln[], int activities[], int birthsteps[]) {return(0)}  /* dummy routine, no effect */
#endif

//------------------------------------------------------- misc ---------------------------------------------------------------------------
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
