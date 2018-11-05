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
#define R_4_enforce2birth 0x10      /* 1: enforce birth for 2 live nbs : choose 1st live nb from top left as ancestor (asym!) */
#define R_5_2uniquepos    0x20      /* 1: for 2-live-nb birth, choose clockwise or anticlockwise (for R_1_choose0nb) elt of canonical bunched pair */
#define R_6_checkmasks    0x40      /* 1: check genetically encoded masks to mask out certain live nb pos's: rule as if these not there */
#define R_7_nongolstat    0x80      /* 1: enforce GoL rule if state of central cell was last changed by a non GoL rule */
#define R_8_nongolstatnbs 0x100     /* 1: enforce GoL rule if state of any cell in nbs was last changed by a non GoL rule */
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
genedata ginitdata = {1,0,-1,0,0,0xfedcba9876543210};  // initialization data structure for gene data
genedata *genedataptr;              // pointer to a genedata instance
HASHTABLE_SIZE_T const* genotypes;  // pointer to stored hash table keys (which are the genotypes)
genedata* geneitems;                // list of genedata structured items stored in hash table
#endif
//-------------------------------------------------------------------------------------------------------------------------------------------------------
int totsteps=0;                     // total number of simulation steps
int statcnts=0;                     // total number of statistic timepts counted
uint64_t codingmask;                // ncoding derived mask for ncoding bits
uint64_t  emptysites = 0;           // cumulative number of empty sites during simulation updates
int nhistG = 0;                     // interval for collecting config histogram data : 0 no collection, nstatG collection with time
int nstatG = 0;                     // interval for collecting other statistical trace data : 0 no collection
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
// get_curgolgstats     get current golgstats array from python
// cmpfunc              compare gene values as numerical unsigned numbers
// cmpfunc1             compare gene counts in population
// countspecies1        count genes with gene array specified as input parameters
// countspecies         count different genes with genes specified at current time point
// cmpfunc2             compare gene values corresponding to given number index in hash table
// cmpfunc3             compare population counts of hash stored genes
// countspecieshash     count different genes in current population from record of all species that have existed
// delay                time delay in ms for graphics
// printxy              terminal screen print of array on xterm
//----------------------------------------------------- begin of subroutines -----------------------------------------------------------------------------

void colorgenes1(uint64_t gol[],uint64_t golg[], uint64_t golgstats[], int cgolg[], int NN2) {
    uint64_t gene, gdiff, g2c, mask;
    int ij,d,d2;
    static int numones[16]={0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};
    if(colorfunction){
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
                        case 10: for (d=0,mask=0;d<16;d++) {d2=(gene>>(d<<2))&0xf;d2=(d2!=0)?1:0;mask+=(d2*0x2a)<<((d%3)<<3);}
                                mask = (mask<<8)+0x7f7f7fff;break;  // 0x2a is 42, approx 1/6 (ie <16/3) of 255
                        case 11: for (d=0,mask=0;d<16;d++) {d2=(gene>>(d<<2))&0xf;d2=numones[d2];mask+=(d2*0xb)<<((d%3)<<3);}
                                mask = (mask<<8)+0x7f7f7fff;break;  // 0x2a is 42, approx 1/6 (ie <16/3) of 255
                        case 12: for (d=0,mask=0;d<64;d++) {d2=(gene>>d)&0x1;mask+=(d2*0x3)<<((d%3)<<3);}
                                mask = (mask<<8)+0x7f7f7fff;break;  // 0x2a is 42, approx 1/6 (ie <16/3) of 255
                        default  : mask = ((d+(d<<6)+(d<<12)+(d<<18))<<8) + 0xff;
                }
                if(colorfunction==2) {
                    if(golgstats[ij]&F_nongolchg) mask = 0x00ffffff;  // color states changed by non GoL rule yellow
                    if(selection>=12) {                               // color as superposition of multiplane gol states
                        for (d=0,mask=0;d<64;d++) {
                            d2=((gol[ij]>>d)&0x1ull);
                            mask+=(d2*0x3)<<((d%3)<<3);
                        }
                        mask = (mask<<8)+0x7f7f7fff;  // 0x2a is 42, approx 1/6 (ie <16/3) of 255
                    }
                    if(selection>=10) {                               // color as superposition of multiplane gol states
                        for (d=0,mask=0;d<16;d++) {
                            d2=((gol[ij]>>(d<<2))&0x1ull);
                            //mask+=(d2*0x2a)<<((d%3)<<3);
                            mask+=(d2*(0x10+((d>>2)<<3)))<<((d%3)<<3);
                        }
                        mask = (mask<<8)+0x7f7f7fff;  // 0x2a is 42, approx 1/6 (ie <16/3) of 255
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
    else{
        // for(d=0;d<256;d++) counts[d]=0;
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
    //for(d=0;d<256;d++)
    //    if (counts[d]) fprintf(stderr,"counting hash table hash %d has %d counts\n",d,counts[d]);
}

void colorgenes(int cgolg[], int NN2) {
    colorgenes1(gol, golg, golgstats, cgolg, NN2);
}

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

extern inline void selectdifft2(uint64_t nb2i, int nb[], uint64_t golg[], uint64_t * birth, uint64_t *newgene) {
// selection based on canonical rotation to bunched pair, choose clockwise or anti-clockwise one (for R_1_choose0nb)
    int k,kmin,nb0,nb1,nbch;
    uint64_t nbmask, nbmaskr, nbmaskrm;
    nb0 = nb2i&0x7;
    nb1 = (nb2i>>4)&0x7;
    nbmask = (0x1ull<<nb0) + (0x1ull<<nb1);
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
        case 0x11ull : k = 4; break;                // 00010001           |...|  |..3|  |0..|   <-
        default  : {                              //                           |...|  |...|
                                                  //                                  |..4|
            fprintf(stderr,"Error in canonical rotation for two live neighbours \nnbmaskrm = %llx\n",nbmaskrm); k = 0;
        } //default case
    } //switch
    if (repscheme & R_1_choose0nb) nbch = nb[kmin];           // replication of live nb in bit 0 of canonical rotation
    else  nbch = nb[(kmin+k)&0x7];                            // replication of live nb in other bit of canonical rotation
    *newgene = golg[nbch];
    *birth = 1ull;
}

extern inline void selectdifft3(uint64_t nbmask, int nb[], uint64_t gol[], int *nbch) {
    int k,kmin;
    uint64_t nbmaskr,nbmaskrm;
    
    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {               // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);                // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                                      // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                        // neighbor mask rotate min is current rotation
            kmin = k;                                                  // no of times rotated to right
        }
    }
    if (repscheme & R_1_choose0nb) *nbch = nb[kmin];                  // replication of live neigbour in bit 0 of canonical rotation
    else {                                                           // replication of live neighbour in most different position
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
        *nbch = nb[(kmin+k)&0x7];                                // rotate unique nb k left (kmin) back to orig nb pat
    }
}

void pack2neighbors(uint64_t gol[], uint64_t golp[]) {           // routine for future use of all 2nd neighbours
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

void update_gol64(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]) {                                           // routine for 64x-gol packed update, not yet used
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

void update_gol16(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]) {
// genes specify on which of 16 planes the gol states at that site are visible: specified by 16 bits 0,4,8,... in gene
// when multiple copy processes are active for gol states, the one with the lowest bit position is applied to most difft gene copying
// gol states looked up for neighbors are the OR of the current plane gol value with gol values on planes specified by the gene
// optionally genes also specify a unique copy plane for a gene: it can only be copied by a rule active on this plane
// in this optional case, copy plane could be lowest non zero bit in gene at pos m for which m%4 == 1 (if none then gene is not copied)
#define deltaxy(ij,x,y)  (ij - (ij&Nmask) + ((ij+(x)&Nmask) + (y)*N)) & N2mask
#ifdef HASH
    genedata gdata;
#endif
    unsigned int ij,ij1,k,kmin,p,p1,nmut,debcnt,mask;
    uint64_t s,su,s3,sg3,s2or3,nbmask,nbmaskr,nbmaskrm,ancestor,newgene,golsh;
    //uint64_t s0,sm;
    uint64_t randnr, randnr2, rand2, statflag;
    int nb1x[8] = {-1,0,1,1,1,0,-1,-1};
    int nb1y[8] = {-1,-1,-1,0,1,1,1,0};
    const uint64_t r1= 0x1111111111111111;
    const uint64_t rc= 0xcccccccccccccccc;

    totsteps++;
    if(!(totsteps%10)) fprintf(stderr,"iteration step %d\r",totsteps);
    if(rulemod) {
      if (selection == 10) {                                                 // gene codes for one secondary plane for each plane  (OLD 10)
        for (ij=0; ij<N2; ij++) {                                            // loop over all sites of 2D torus with side length N
            for (p=0,golsh=0ull;p<16;p++) {                                  // 16 different planes, one every 4 bits in gol,golg
                p1 = (golg[ij]>>(p<<2))&0xf;
                golsh |= ((gol[ij]>>(p1<<2))&1ull)<<(p<<2);                  // collect active bits from 2nd planes pointed to by gene
                // golsh |= (gol[ij]>>((golg[ij]>>(p<<2))<<2))<<(p<<2);      // erroneous earlier version of commit 61a7463 giving exploration waves
            }
            golmix[ij] = golsh/*|gol[ij]*/;                                  // bits active from secondary plane (or primary if or with gol)
        }
      }
      else { // selection == 11                                              // gene codes for 4(2) nearest plane neighbour masks
        for (ij=0; ij<N2; ij++) {                                            // loop over all sites of 2D torus with side length N
            for (p=0,golsh=0ull;p<16;p++) {                                  // 16 different planes, one every 4 bits in gol,golg
                mask = (golg[ij]>>(p<<2))&0xf;                               // bit mask indicating which nearby planes visible from plane
                for (k=1;k<3;k++) {                                          // 2 nearest planes variant
                //for (k=0;k<4;k++) {                                        // 4 nearest planes variant
                    if((mask>>k)&0x1) {
                        p1 = (p+k+(k>>1)-2)&0xf;
                        golsh ^= ((gol[ij]>>(p1<<2))&0x1ull)<<(p<<2);        // nearby planes p1=p+{-2,-1,1,2} for k 0,1,2,3   or p+{-1,1} for k 1,2
                    }
                }
            }
            golmix[ij] = golsh/*|gol[ij]*/;                                  // bits active from (primary or) two(four) secondary planes
        }
      }
    }
    else { // no gol rule modification, independent planes
         for (ij=0; ij<N2; ij++) golmix[ij] = gol[ij];                       // loop over all sites of 2D torus with side length N
    }
    // for (ij=0; ij<N2; ij++) if ((golmix[ij]&r1) != golmix[ij]) fprintf(stderr,"error in golmix calculation %llx",golmix[ij]);
    
    for (ij=0; ij<N2; ij++) {                                                   // loop over all sites of 2D torus with side length N
        for(k=0,s=0ull;k<8;k++) {                                               // compute no of live neighbours using golmix (coupled plane neighbours)
            ij1=deltaxy(ij,nb1x[k],nb1y[k]);
            s +=golmix[ij1];
        }
        su = s&rc;                                                              // upper 2 bits (3,2) of sum of live neighbours : non zero if sum 4-8
        sg3 = (((su>>1)|su)>>2) & r1;                                           // for each plane 1 if sum is gt 3
        s2or3 = (~sg3)&(s>>1)&r1;                                               // sum is 2 or 3 for each plane, ie not gt 3 and bit 2 is 1
        s3 = s2or3&s&r1;                                                        // sum is 3 for each plane (s2or3 and bit 0 is 1)
        newgol[ij]=s2or3&(gol[ij]|s3)&r1;                                       // parallel gol rule with coupled sum s for each plane
        statflag = 0ull;
        if((golmix[ij]^gol[ij])&s2or3) statflag |= F_notgolrul;                 // not gol rule in at least one plane
        if (s3) {                                                               // birth in at least one plane
          for(p=0;p<16;p++) {
            if((s3>>(p<<2))&1ull) {
                for(k=0,nbmask=0ull;k<8;k++) {
                    ij1=deltaxy(ij,nb1x[k],nb1y[k]);
                    nbmask |= (((golmix[ij1])>>(p<<2))&1ull)<<k;
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
                nmut = (randnr >> 56) & 0x3f;                                     // choose mutation position for length 64 gene : from bits 56:61 of randnr
                ancestor = newgene;
                newgene = newgene ^ (rand2<<nmut);                                // introduce single mutation with probability pmut = probmut
                newgolg[ij]=newgene;
                statflag = statflag | F_birth;
                if (rand2) statflag = statflag | F_mutation;
                newgolgstats[ij] = statflag;
#ifdef HASH
                if((genedataptr = (genedata *) hashtable_find(&genetable, golg[ij])) != NULL) { // central old gene present: overwritten
                    genedataptr->popcount--;  // if 0 lastextinctionframe updated after whole frame calculated
                }
                else fprintf(stderr,"step %d hash storage error 1 in update_gol16, gene %llx at %d not stored\n",totsteps,golg[ij],ij);
                if((genedataptr = (genedata *) hashtable_find(&genetable, newgene)) != NULL) {
                    genedataptr->popcount++;
                }
                else {
                    gdata=ginitdata;
                    gdata.firstbirthframe = totsteps;
                    gdata.firstancestor = ancestor;
                    hashtable_insert(&genetable, newgene,(genedata *) &gdata);
                }
#endif
                p=16;                                                             // break from loop at first birth
            } //if 3 live nbs at p
        } // for p
      }  // if (s3)
      else {
        newgolg[ij]=gene0;
#ifdef HASH
        if(golg[ij]!= gene0) {
            newgene= gene0;                                                     // default is relax to uncoupled gene

            if((genedataptr = (genedata *) hashtable_find(&genetable, golg[ij])) != NULL) { // central old gene present: overwritten
                    genedataptr->popcount--;  // if 0 lastextinctionframe updated after whole frame calculated
            }
            else fprintf(stderr,"step %d hash storage error 1 in update_gol16, gene %llx at %d not stored\n",totsteps,golg[ij],ij);
            if((genedataptr = (genedata *) hashtable_find(&genetable, newgene)) != NULL) {
                genedataptr->popcount++;
            }
            else {
                gdata=ginitdata;
                gdata.firstbirthframe = totsteps;
                gdata.firstancestor = gene0;
                hashtable_insert(&genetable, newgene,(genedata *) &gdata);
            }
        }
#endif
      } // end else (s3)
    } // for ij
    for (ij=0; ij<N2; ij++) {        // update lattices
        gol[ij] = newgol[ij];        // copy new gol config to old one
        golg[ij] = newgolg[ij];      // copy new genes to old genes
        golgstats[ij] = newgolgstats[ij]; // copy new golg stat flags to old genes
    }
#ifdef HASH
    for (ij=0; ij<N2; ij++) {       // complete missing hash table records including activities
        if((genedataptr = (genedata *) hashtable_find(&genetable, golg[ij])) != NULL) {
            if(genedataptr->popcount == 0) {
                genedataptr->lastextinctionframe = totsteps;
                genedataptr->nextinctions++;
            }
        }
        else fprintf(stderr,"step %d hash storage error 2 in update_gol16, gene %llx at %d not stored\n",totsteps,golg[ij],ij);
        if((genedataptr = (genedata *) hashtable_find(&genetable, newgolg[ij])) != NULL)
            genedataptr->activity ++;
        else fprintf(stderr,"step %d hash storage error 3 in update_gol16, gene %llx at %d not stored\n",totsteps,golg[ij],ij);
    }
#endif

}

void update(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]){
    /* update GoL for toroidal field which has side length which is a binary power of 2 */
    /* encode without if structures for optimal vector treatment */

    int s, s0, k, k1, nmut;
    unsigned int checkmasks,checknongol,rulemod1,rulemod2,rulemod1ij,rulemod2ij;
    int nb[8], nbc, nbch, ij, i, j, jp1, jm1, ip1, im1;
    uint64_t g, gs, nb1i, nb2i, randnr, randnr2, r2;
    uint64_t nbmask, nbmaskr;
    uint64_t newgene, ancestor, livegenes[3];
    uint64_t s2or3, birth, statflag, nextgolstate;
#ifdef HASH
    genedata gdata;
#endif

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
                  selectdifft3(nbmask, nb, gol, &nbch);
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
                        else if (repscheme & R_5_2uniquepos) selectdifft2(nb1i,nb,golg,&birth,&newgene);
                        else selectone(s,nb1i,nb,golg,&birth,&newgene);
                        if(!birth && (repscheme&R_4_enforce2birth)) selectdifft2(nb1i,nb,golg,&birth,&newgene);
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
#ifdef HASH
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
#endif
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
#ifdef HASH
                    if(gol[ij]) {                                           // death : need to update hash table
                        if((genedataptr = (genedata *) hashtable_find(&genetable, golg[ij])) != NULL) {
                            genedataptr->popcount--;
                        }
                        else fprintf(stderr,"hash storage error 2, gene %llx not stored\n",golg[ij]);
                    }
#endif
                    newgol[ij]  = 0ull;                                     // new game of life cell value dead
                    newgolg[ij] = 0ull;                                     // gene dies or stays dead
                    if(gol[ij]) statflag |= F_death;
                }
            } // end no birth
        }  // end if s2or3
        else {                                                              // else not birth or survival, 0 values for gol and gene
#ifdef HASH
            if(gol[ij]) {                                                   // death : need to update hash table
                if((genedataptr = (genedata *) hashtable_find(&genetable, golg[ij])) != NULL) {
                    genedataptr->popcount--;
                }
                else fprintf(stderr,"hash storage error 3, gene %llx not stored\n",golg[ij]);
            }
#endif
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
        emptysites = emptysites + newgol[ij];                               // keep track of empty sites, same information as total activity of occupied sites
        newgolgstats[ij] = statflag;
    }  // end for ij

#ifdef HASH
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
#endif

    for (ij=0; ij<N2; ij++) {        // update lattices
	    gol[ij] = newgol[ij];        // copy new gol config to old one
	    golg[ij] = newgolg[ij];      // copy new genes to old genes
        golgstats[ij] = newgolgstats[ij]; // copy new golg stat flags to old genes
    }

}

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

void genelife_update (int nsteps, int nhist, int nstat) {
    /* update GoL for toroidal field which has side length which is a binary power of 2 */
    /* encode without if structures for optimal vector treatment */
    int t;
    uint64_t *newgol, *newgolg;

    nhistG = nhist;
    nstatG = nstat;
    for (t=0; t<nsteps; t++) {
        newgol = planes[newPlane];
        newgolg = planesg[newPlane];

        if (selection<10) update(gol,golg,newgol,newgolg);                    // calculate next iteration with selection
        else if (selection <12) update_gol16(gol,golg,newgol,newgolg);        // calculate next iteration for 16x multiplane version
        else update_gol64(gol,golg,newgol,newgolg);                           // calculate next iteration for 64x multiplane version
        if(nhist && (totsteps%nhist == 0)) countconfigs();                    // count configurations
        if(nstat && (totsteps%nstat == 0)) tracestats(gol,golg,golgstats,N2); // time trace point

        curPlane = (curPlane +1) % numPlane;            // update plane pointers to next cyclic position
        newPlane = (newPlane +1) % numPlane;
        gol = planes[curPlane];                         // get planes of gol,golg data
        golg = planesg[curPlane];
    }
}

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
    int ij,ij1,i0,j0,i,j,Nf,k,cnt,icf,p,nstartgenes;
    uint64_t g;

#ifdef HASH
    int hcnt;
#endif

    static unsigned int rmask = (1 << 15) - 1;
    static int notfirst = 0;
    // Range: rand returns numbers in the range of [0, RAND_MAX ), and RAND_MAX is specified with a minimum value of 32,767. i.e. 15 bit

    uint64_t startgenes[16];
    char *golgin;
    
    srand(1234567);
    state[0] = rand();state[1] = rand();
    cnt = 0;
    totsteps = 0;
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
    
    if (selection == 10) gene0=0xfedcba9876543210;
    else gene0=0x0ull;
    
    nstartgenes = 8;
    if (selection >= 10) nstartgenes = 16;
    switch (selection) {
        case 0: for (k=0;k<4;k++) {startgenes[k]=0xf0f0f0f0f0f0f0f0;startgenes[k+4]=0x0f0f0f0f0f0f0f0f;} break;
        case 1: for (k=0;k<8;k++) startgenes[k]=((0x1ull<<k*3)-1ull)<<20;break;
        case 2:
        case 3: for (k=0;k<8;k++) startgenes[k]=(((0x1ull<<20)-1ull)<<20)+((0x1ull<<k)-0x1ull);break;
        case 4: for (k=0;k<8;k++) {g = 0xfffff0ull + k; startgenes[k] = k<4 ? g : ~g;} break;
        case 5: for (k=0;k<8;k++) {g = 0xf0ull + k; startgenes[k]= k<4 ? g : (~g)|(0xfull<<16);} break;
        case 10: for (p=0,g=0ull;p<16;p++) g|= (p&0xfull)<<(p<<2);for (k=0;k<16;k++) startgenes[k] = g; // first set all startgenes as uncoupled
         for (k=0;k<16;k++) startgenes[k] = (startgenes[k] & (~(0xfull<<(k<<2)))) | ((k+1)&0xfull)<<(k<<2);break; // couple plane k+1 to k in startgene k
        case 11: for (k=0;k<16;k++) {for (p=0,g=0ull;p<16;p++) g|= (k&0xfull)<<(p<<2);startgenes[k] = g;} break;
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
            golg[ij] = 0;
            if (selection == 10) golg[ij]=0xfedcba9876543210;
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
                gol[ij] = 0;
                golg[ij] = 0;
                if (selection == 10) golg[ij]=gene0;
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
            golg[ij] = 0ull;
            if (selection == 10) golg[ij]=gene0;
            golgstats[ij] = 0;
        }
        i0 = j0 = (N>>1)-(Nf>>1);
        for (i=0; i<Nf; i++) {
            for (j=0; j<Nf; j++) {
                ij=i0+i+N*(j0+j);
                if(selection<10) gol[ij] = ((rand() & rmask) < initial1density)?1ull:0ull;
                else for (k=0;k<16;k++) gol[ij] |= ((rand() & rmask) < initial1density)?(1ull<<(k<<2)):0ull;
            }
        }
        for (ij=0; ij<N2; ij++) {
            g = 0ull;
            if (gol[ij]||(selection>=10)) { // if live cell or multiplane, fill with random genome g or randomly chosen startgene depending on initialrdensity
                if (((unsigned) rand() & rmask) < initialrdensity) for (k=0; k<64; k++) g = (g << 1) | (rand() & 0x1);
                else if (startgenechoice == nstartgenes) g = startgenes[0xf & rand() & (nstartgenes-1)];
                else if (startgenechoice > nstartgenes) fprintf(stderr,"startgenechoice %d out of range\n",startgenechoice);
                else g = startgenes[0xf & startgenechoice & (nstartgenes-1)];
                cnt++;
            }
            golg[ij] = g;
            // if (golg[ij] == 0ull && gol[ij] != 0ull) fprintf(stderr,"zero gene at %d\n",ij);
        }
        // for (ij=0; ij<40; ij++) fprintf(stderr,"gene at %d %llx\n",ij,golg[ij]);   // test first 40
    }

#ifdef HASH
    for (ij=0; ij<N2; ij++) {
        if(gol[ij]||(selection>=10)) {
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
#endif
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

void get_curgolgstats(uint64_t outgolgstats[], int NN){
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolgstats[ij] = golgstats[ij];                       // Note that golgstats is not dealt with in planes !
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
    fprintf(stderr,"at step %d cumulative activity = %llu\n",totsteps,(N2 * (uint64_t) totsteps) - emptysites);
    fprintf(stderr,"rulemod\trepscheme\tselection\toverwritemask\tsurvival\n");
    fprintf(stderr,"%d\t%d\t\t%d\t\t%d\t\t%d\n",rulemod,repscheme,selection,overwritemask,survival);
    fprintf(stderr,"nlog2pmut\tinit1\tinitr\tncoding\tstartchoice\n");
    fprintf(stderr,"%d\t\t%d\t%d\t%d\t%d\n",nlog2pmut,initial1density,initialrdensity,ncoding,startgenechoice);
}

void countspecies() {
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
    // fprintf(stderr,"cumulative activity = %llu\n",(N2 * (uint64_t) totsteps) - emptysites);
}
#else
void countspecieshash() {}  /* dummy routine, no effect */
#endif

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
