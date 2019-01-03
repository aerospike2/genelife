// subgenelife.c
// Written by John S. McCaskill and Norman H. Packard
//
// Project fastgenegol
//
// First created by John McCaskill on 14.07.17. Last modified Nov 2018.
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
unsigned int rulemod = 1;           // det: whether to modify GoL rules
int selection = 1;                  // fitness of 2 live neighbours:
                                    // 0. integer value   1. number of ones  2. scissors-stone-well-paper: wins over left 1-1-2-2
                                    // 3. scissors-stone-well-paper 1-1-1-1 4. 5. predator prey 6. two target coding 7. selection always fails (no result)
                                    // 8. genes encode lut based on sum for survival (1-8) and birth (1-8)
                                    // 9. like 8 but genes penalized in fitness for the number of lut entries active (& poss. also proximity to mask edge)
                                    // 10. genes encode lut based on sum and canonical rotations survival/birth sum 2-6 (4,7,10,7,4 canon. rotations)
                                    // 11. like 10 but with genes penalized as in 9
                                    // 12. genes encode lut based on total live neighbour distance lower and upper cutoffs for survival and birth
                                    // 13. like 12 but with genes penalized as in 9
                                    // 14. two match-coupled planes : one with genes
                                    // 15. like 14 but with genes penalized as in 9
                                    // 16. 1-16 planes coupled to genes via plane index
                                    // 18. 1-16 planes coupled to 2-4 nearest planes via gene masks
                                    // 20. 64 planes of GoL currently independent
unsigned int repscheme = 1;         // replication scheme: lowest 10 bits define 5 pairs of bit options for:
                                    // 0,1 select birth 2,3 neighbour choice 4,5 enforce birth 6,7 2nd,1st neighbour genes 8,9 no successive nonGoL
                                    // bits 10-13 specify 4-bit mask for 2-live neighbour birth onl for those of the 4 canonical configurations set
                                    // bits 14-20 activate quadrant exploration of specific subset of the first 5 pairs of repscheme and survival and overwrite masks
unsigned int survivalmask = 0x3;    // for selection 0-7 survive mask for two (bit 1) and three (bit 0) live neighbours
                                    // for selection=8,9 it is 16-bit and for 10-11 it is 19-bit and 12-13 32-bit survival mask to restrict luts
unsigned int birthmask = 0x0;       // for selection 0-7 unused, birth control is done via bits 0,1,4,5 of repscheme
                                    // for selection=8,9 it is 16-bit and for 10-11 it is 19-bit and 12-13 32-bit birthmask to restrict luts
unsigned int overwritemask = 0x3;   // for selection 0-7 bit mask for 4 cases of overwrite: bit 0. s==3  bit 1. special birth s==2
                                    // for selection=8,9 it is 16-bit and for 10-11 it is 19-bit and 12-13 32-bit birthmask to restrict luts

int ncoding = 1;                    // byte 0 of python ncoding : number of coding bits per gene function
int ncoding2 = 0;                   // byte 1 of python ncoding: number of coding bits per gene function for masks in connection with repscheme add2ndmask1st R_6,7
int NbP;                            // byte 2 of python ncoding: number of bit planes used for parallel gol planes packed into gol long unsigned integer : used for selection = 14+
int NbG;                            // number of bits in gene used : fully determined by selection and NbP
unsigned int pmutmask;              // binary mask so that prob of choosing zero is pmut = pmutmask/2^32. If value<32 interpret as integer -log2(prob).
//-----------------------------------------------------------initialization and color parameters---------------------------------------------------------
int initial1density = (1<<15)>>1;   // initial density of ones in gol as integer value, divide by 2^15 for true density
int initialrdensity = (1<<15)>>1;   // initial density of random genes in live sites, divide by 2^15 for true density
int startgenechoice = 8;            // selection for defined starting genes 0-8 (8 is random 0-7) otherwise choose only particular startgene
int initfield = 0;                  // 0 input from random field of random or start genes, 1 input from file genepat.dat of 32x32 indexes to start genes
                                    // value i>1: initialized to random field on central ixi block only, outside this zero.
int colorfunction = 0;              // color function choice of 0: hash or 1: functional (color classes depends on selection parameter)
                                    // 2: as in 1 but color sites where last step was non GoL rule yellow, 3: as in 2 but yellow if state produced by non GoL
                                    // 4: activities 5: genealogies without time 6: genealogies with time 7: genealogies with time and activity 8: gliders
#define ASCII_ESC 27                /* escape for printing terminal commands, such as cursor repositioning : only used in non-graphic version */
//-----------------------------------------masks for named repscheme bits (selection 0-7) ----------------------------------------------------------------
#define R_0_2sel_3live     0x1      /* 1: for 3-live-n birth, employ selection on two least different live neighbours for ancestor */
#define R_1_2sel_2live     0x2      /* 1: allow 2-live-n birth, employ selection on 2 live neighbours for ancestor */
#define R_2_canonical_nb   0x4      /* 1: choose live neighbour at zero bit in canonical rotation 0: choose most difft position */
#define R_3_neutral_pos    0x8      /* 1: for 2-live-nb birth, choose canonical position as specified by R_2_canonical_nb rather than doing selection */
#define R_4_enforce3birth  0x10     /* 1: enforce birth for 3 live nbs (with ancestor most difft) in case 2-select fails (req. R_0_2sel_3live==1 for effect) */
#define R_5_enforce2birth  0x20     /* 1: enforce birth for 2 live nbs : choose 1st live nb from top left of canonical config as ancestor (asym!) */
#define R_6_2ndnb_genes    0x40     /* 1: execute genetically encoded 2nd neighbours of live neighbours to determine birth & 1st nb ancestor */
#define R_7_1stnb_masks    0x80     /* 1: check genetically encoded masks to mask out certain live nb pos's: rule as if these not there */
#define R_8_nongolstat     0x100    /* 1: enforce GoL rule if state of central cell was last changed by a non GoL rule */
#define R_9_nongolstatnbs  0x200    /* 1: enforce GoL rule if state of any cell in nbs was last changed by a non GoL rule */
#define R_10_2birth_k0     0x400    /* 1: bit position at start of k1-4 mask for selective subset of 2-births */
#define R_10_13_2birth_k4  0x3c00   /* 1: enforce birth for 2 live nbs canonical config for one of k= 1,2,3,4, next 4 bits: choose 1st live nb from TL (asym!) */
#define R_quadrant         0x1fc000 /* 1: quarter the spatial domain with one or more of 7 pairs of repscheme bits ie 4 different values */
#define R_14_quadrant_sele 0x4000   /* q0 1: quarter the spatial domain with selection enable values for 2,3 live nbs: only in update ie for selection<8 */
#define R_15_quadrant_posn 0x8000   /* q1 1: quarter the spatial domain with selection enable values for 2,3 live nbs: only in update ie for selection<8 */
#define R_16_quadrant_enfb 0x10000  /* q2 1: quarter the spatial domain with enforce birth values for 2,3 live nbs: only in update ie for selection<8 */
#define R_17_quadrant_2nb1 0x20000  /* q3 1: quarter the spatial domain with 1st nb masks and/or 2nd nb addition: only in update ie for selection<8 */
#define R_18_quadrant_ngol 0x40000  /* q4 1: quarter the spatial domain with last non gol rule and/or non gol created state: only in update ie for selection<8 */
#define R_19_quadrant_surv 0x80000  /* q5 1: quarter the spatial domain with survival values for 2,3 live nbs: only in update ie for selection<8 */
#define R_20_quadrant_over 0x100000 /* q6 1: quarter the spatial domain with overwrite values for 2,3 live nbs: only in update ie for selection<8 */
//.................................................. LUT repscheme (selection 8-15) repscheme bits .......................................................
#define R_0_survivalgene  0x1       /* 1: survival gene chosen from central existing gene 0: survival gene taken from neighbours as in birth */
#define R_1_nb_OR_AND     0x2       /* 1: OR of neighbours determines genetic LUT in selection 8,10,12,14 0: AND of neighbours */
#define R_2_canonical_nb  0x4       /* 1: choose live neighbour at zero bit in canonical rotation 0: choose most difft position */
#define R_46_repselect    0x70      /* 0-7 choice of selection mechanism for LUT genes : 0: min 1: max 2: min 1s 3: max 1s 4: neutral 5: neutral difft 6,7: c-S-2B */
#define R_7_random_resln  0x80      /* 1: random choice amongst selected live neighbours 0: deterministic choice based on gene content and position */
#define R_810_disambig    0x700     /* 0-7 choice of different disambiguation mechanisms for symmetric canonical rotations */
//................................................... multiplane (selection 16-19) repscheme bits ........................................................
#define RP_0_or_current   0x1       /* 1: do or between coupling plane and current plane(s) 0: just take coupling plane(s) */
#define RP_1_c4planes     0x2       /* 1: 4 nearest planes coupled 0: 2 nearest planes coupled */
#define RP_2_survive3     0x4       /* 1: survival for 3 live-nbs 0: no survival for 3 live-nbs */
#define RP_3_survive2     0x8       /* 1: survival for 2 live-nbs 0: no survival for 2 live-nbs */
#define RP_4_grow_planes  0x10      /* 1: grow from plane 0 0: chose coupling to any plane */
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
//.......................................................................................................................................................
hashtable_t quadtable;              // hash tabel for quad tree
int quadtree = 0;                   // whether to do quadtree calculation
struct quadkey;
typedef struct quadkey {
    struct quadkey *next;
    uint64_t key;
} quadkey;
struct quadnode;                    // prototype for quadtree node defined below with recursive elements
typedef union quad {                // node pointer and 64-bit patterns occupy the same memory : access Union elts as quad.node or quad.patt
    struct quadnode *node;
    uint64_t patt;
} quad;
typedef struct quadnode {           // stored quadtree binary pattern nodes for population over time (currently only for analysis not computation)
   struct quadkey *next;            // hash-linked list of elements starting with same key : follow hash key chain
   quad nw, ne, sw, se;             // constant pointers to quadnodes or 64-bit patterns : we terminate one level higher than Gosper & golly
   struct quadnode *res;            // store cached value associated with the node (optional here)
   unsigned char isnode,active;     // 1 if this is a node not a pattern and an active part of the current time step CA resp.
   unsigned short reserve;
   unsigned int hits;               // number of references to finding this pattern
   unsigned int pop1s;              // 32 bit number of 1s in quadnode
   unsigned int firsttime;          // first time node is identified
} quadnode;
quadnode quadinit = {NULL,0ull,0ull,0ull,0ull,NULL,0,1,0,1,0,0};
quadnode * qimage;
quadkey *freenodes;
int quadchainmem = 0;
int quadcollisions = 0;
int nquadpatterns = 0;
#endif                              // HASH
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
uint64_t *golgstats, *newgolgstats; // pointers to 64 bit masks for different events during processing at one of the plane cycle locations
uint64_t golmix[N2];                // array for calculating coupling of genes between planes for multiplane genelife, or for packing configs
uint64_t gene0;                     // uncoupled planes background gene, non zero for selection==16,17
uint64_t selectedgene;              // gene currently selected interactively in graphics window
uint64_t genegol[16];               // genes encoding for GoL in various LUT selection models indexed as selection-8
//---------------------------------------------------------additional variables set in simulation--------------------------------------------------------
unsigned int canonical;             // current value of choice of canonical position repscheme bit 2 : needed globally in ...difft2-6 routines
int quadrants=-1;                   // integer choice of bit pair from repscheme/survivalmask/overwritemask for quadrant division of array (-1 none)
uint64_t displayplanes;             // mask for displaying planes in selection models 16+
uint64_t displayoneplane;           // display one plane only if <64, all planes if 64
int randomsoup=0;                   // 1,2 steady or intermittent random input of genes and gol states into the square central region defined by initfield
int rbackground=0;                  // integer background rate of random gene input per frame per site : rate is rbackground/32768 (nonzero overides soup)
                                    // the gene input depends on randomsoup value: 2 GoL 1 random
int vscrolling=0;                   // whether to do vertical scrolling to track upwards growth (losing all states that fall off downward cliff)
int last_scrolled = 0;              // whether vscrolling applied on last time step (needed for correct glider detection)
int ymax = 2000;                    // activity scale max for plotting : will be adjusted dynamically or by keys
double log2ymax = 25.0;             // activity scale max 2^25 = 33.5 * 10^6
int activitymax;                    // max of activity in genealogical record of current population
//------------------------------------------------ arrays for time tracing, activity and genealogies ----------------------------------------------------
const int startarraysize = 1024;    // starting array size (used when initializing second run)
int arraysize = startarraysize;     // size of trace array (grows dynamically)
int *livesites = NULL;              // dynamic array pointer for statistics of number of live sites over time
int *genestats = NULL;              // dynamic array pointer for statistics of number of 4 genotype classes over time
int *stepstats = NULL;              // dynamic array pointer for statistics of site update types over time
int *configstats = NULL;            // dynamic array pointer for statistics of gol site configurations (x,y,t) offsets
uint64_t acttrace[N2];              // scrolled trace of last N time points of activity gene trace
uint64_t genealogytrace[N2];        // image trace of genealogies for N most frequently populated genes
uint64_t working[N2];               // working space array for calculating genealogies and doing neighbour bit packing
//.......................................................................................................................................................
short unsigned int label[N2];       // labels for connected component labelling
typedef struct equivrec {           // equivalence record
  short unsigned int parent;
  short unsigned int rank;
  // short unsigned int size;
} equivrec;
equivrec equiv[N2<<2];              // equivalences between labels
typedef struct component {
    short unsigned int N,S,E,W;
    unsigned int first;
    unsigned int corner;
    unsigned int Nside;
    quadnode * quad;
} component;
component complist[N2<<2];
int ncomponents = 0;
//------------------------------------------------ planes and configuration offsets----------------------------------------------------------------------
int offdx=0,offdy=0,offdt=0;        // display chosen offsets for glider analysis with colorfunction 8
int Noff = 9;                       // number of offsets
int **offsets;                      // array of offsets (2D + time) for planes
int *histo;
int numHisto;
// initialize planes:
#define maxPlane 4                  /* maximum number of planes allowed : values 2,4,8 allowed */
int curPlane = 0;                   // current plane index
int newPlane = 1;                   // new plane index
int numPlane = maxPlane;            // number of planes must be power of 2 to allow efficient modulo plane
uint64_t *planes[maxPlane];         // ring buffer planes of gol array states
uint64_t *planesg[maxPlane];        // ring buffer planes of golg genes
uint64_t *planesgs[maxPlane];       // ring buffer planes of golgstatus bits
uint64_t plane0[N2];                // gol   0
uint64_t plane1[N2];                // gol   1
uint64_t planeg0[N2];               // golg  0
uint64_t planeg1[N2];               // golg  1
uint64_t planegs0[N2];              // golgs 0
uint64_t planegs1[N2];              // golgs 1
#if maxPlane >= 2
uint64_t plane2[N2];                // gol   2
uint64_t plane3[N2];                // gol   3
uint64_t planeg2[N2];               // golg  2
uint64_t planeg3[N2];               // golg  3
uint64_t planegs2[N2];              // golgs 2
uint64_t planegs3[N2];              // golgs 3
#endif
#if maxPlane >= 4
uint64_t plane4[N2];                // gol   4
uint64_t plane5[N2];                // gol   5
uint64_t plane6[N2];                // gol   6
uint64_t plane7[N2];                // gol   7
uint64_t planeg4[N2];               // golg  4
uint64_t planeg5[N2];               // golg  5
uint64_t planeg6[N2];               // golg  6
uint64_t planeg7[N2];               // golg  7
uint64_t planegs4[N2];              // golgs 4
uint64_t planegs5[N2];              // golgs 5
uint64_t planegs6[N2];              // golgs 6
uint64_t planegs7[N2];              // golgs 7
#endif
//------------------------------------------------------- fast macro random number generator ------------------------------------------------------------
                                    // Wikipedia "Xorshift" rewritten here as inline macro &
                                    // Vigna, Sebastiano. "xorshift*/xorshift+ generators and the PRNG shootout". Retrieved 2014-10-25.
static uint64_t state[2];           // State for xorshift pseudorandom number generation. The state must be seeded so that it is not zero
#define RAND128P(val) {                                                       \
    uint64_t x = state[0]; uint64_t const y = state[1];                       \
	state[0] = y;	x ^= x << 23;  state[1] = x ^ y ^ (x >> 17) ^ (y >> 26);  \
	val = state[1] + y;}
//.......................................................................................................................................................
const uint64_t m1  = 0x5555555555555555; //binary: 0101...           Constants for Hamming distance macro POPCOUNT24C
const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...
#define POPCOUNT64C(x, val) {                  /* Wikipedia "Hamming Weight" popcount4c alg */  \
    uint64_t xxxx;                             /* define copy of x argument so that we do not change it */ \
    xxxx = x;                                  /* copy x argument */ \
    xxxx -= (xxxx >> 1) & m1;                  /* put count of each 2 bits into those 2 bits */ \
    xxxx = (xxxx & m2) + ((xxxx >> 2) & m2);   /* put count of each 4 bits into those 4 bits */ \
    xxxx = (xxxx + (xxxx >> 4)) & m4;          /* put count of each 8 bits into those 8 bits */ \
    val = (xxxx * h01) >> 56;}                 /* left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... */
//.......................................................................................................................................................
#define PATTERN2_32(x, pat, found) {           /* find number of 2-bit aligned 2-bit copies of pattern pat in 32-bits of 64-bit integer x */ \
    const uint64_t r1_32 = 0x11111111ull;       \
    const uint64_t r3_32 = 0x33333333ull;       \
    const uint64_t r5_32 = 0x55555555ull;       \
    const uint64_t rf_32 = 0xffffffffull;       \
    uint64_t xxxx;                             /* define working variable */ \
    xxxx=(uint64_t) pat;                       /* copy pat to ensure it is not assumed to be a variable that can be changed */ \
    xxxx|=xxxx<<2;                             /* doubles the pattern to the left */ \
    xxxx|=xxxx<<4;                             /* doubles the patterns to the left : now 4 copies */ \
    xxxx|=xxxx<<8;                             /* doubles the patterns to the left : now 8 copies */ \
    xxxx|=xxxx<<16;                            /* doubles the patterns to the left : now 16 copies */ \
    xxxx = (xxxx ^ (uint64_t) x) & rf_32;      /* xor x argument with 16 copies of pat */ \
    if (xxxx) {                                /* need to exclude the 16 match case as the fast count fails here */ \
        xxxx = ~xxxx&rf_32;                    /* invert difference map to yield identity map, 2 ones if match at a position */ \
        xxxx = (xxxx&(xxxx>>1))&r5_32;    /* convert 2-bit set of identity patterns to 1/0 decision at bit 0 of 2 bit pattern */ \
        xxxx = (xxxx+(xxxx>>2))&r3_32;    /* do first sum pairwise so that set of 8 2-bit sums spaced by 4-bits */ \
        found = ((xxxx&r1_32) * r1_32) >> 28;} /* found is returned as the number of patterns found at any of the 16 positions : <16 */ \
    else found = 16;}                          /* match at all positions */
//.......................................................................................................................................................
const uint64_t r1 = 0x1111111111111111ull;
#define PATTERN4(x, pat, found) {              /* find number of 4-bit aligned 4-bit copies of pattern pat in 64-bit x */ \
    uint64_t xxxx;                             /* define working variable */ \
    xxxx=(uint64_t) pat;                       /* copy pat to ensure it is not assumed to be a variable that can be changed */ \
    xxxx|=xxxx<<4;                             /* doubles the pattern to the left */ \
    xxxx|=xxxx<<8;                             /* doubles the patterns to the left : now 4 copies */ \
    xxxx|=xxxx<<16;                            /* doubles the patterns to the left : now 8 copies */ \
    xxxx|=xxxx<<32;                            /* doubles the patterns to the left : now 16 copies */ \
    xxxx ^= (uint64_t) x;                      /* xor x argument with 16 copies of pat */ \
    if (xxxx) {                                /* need to exclude the 16 match case as the fast count fails here */ \
        xxxx = ~xxxx;                          /* invert difference map to yield identity map, 4 ones if match at a position */ \
        xxxx &= (xxxx>>1)&(xxxx>>2)&(xxxx>>3); /* convert 4-bit set of identity patterns to 1/0 decision at bit 0 of 4 bit pattern */ \
        found = ((xxxx&r1) * r1) >> 60;}       /* found is returned as the number of patterns found at any of the 16 positions : <16 */ \
    else found = 16;}                          /* match at all positions */
//.......................................................................................................................................................
// const uint64_t h01 = 0x0101010101010101;    /* multiplicand defined already for POPCOUNT24: used to sum up all 8 8-bit bytes */
#define PATTERN8(x, pat, found) {              /* find number of 8-bit aligned 8-bit copies of pattern pat in 64-bit x */ \
    uint64_t xxxx;                             /* define working variable */ \
    xxxx=(uint64_t) pat;                       /* copy pat to ensure it is not assumed to be a variable that can be changed */ \
    xxxx|=xxxx<<8;                             /* doubles the pattern to the left : now 2 copies */ \
    xxxx|=xxxx<<16;                            /* doubles the patterns to the left : now 4 copies */ \
    xxxx|=xxxx<<32;                            /* doubles the patterns to the left : now 8 copies */ \
    xxxx ^= (uint64_t) x;                      /* xor x argument with 8 copies of pat */ \
    xxxx = ~xxxx;                              /* invert difference map to yield identity map, 8 ones if match at a position */ \
    xxxx &= xxxx>>1;                           /* convert 8-bit set of identity patterns to 1/0 decision */ \
    xxxx &= xxxx>>2;                           /*    at bit 0 of 8 bit pattern */ \
    xxxx &= xxxx>>4;                           /*    in 3 steps */  \
    found = ((xxxx&h01) * h01) >> 56;}         /* found is returned as the number of patterns found at any of the 8 positions */
//.......................................................................................................................................................
#define FIRST1INDEX(v, c) {                    /* starting point 64bit from Sean Eron Anderson https://graphics.stanford.edu/~seander/bithacks.html#ZerosOnRightParallel */  \
    uint64_t mmmm,mmmq;                        /* note also arguments must be of types uint64_t and int respectivley */ \
    int cccc;                                  /* takes on successive integer values 32,16,84,2,1 */ \
    int tttt;                                  /* logical to integer variable true=one false=zero : if a 1 in v under mask mmmm */ \
    c=v?0:1;                                   /* c will contain count of number of zeros on right of last one, here if v is all zeros then start from 1 */ \
    mmmm=~0ull;                                /* all ones, mask from previous stage */ \
    for (cccc=1<<5;cccc>0;cccc>>=1) {          /* loop over cccc goes 32,16,8,4,2,1 */ \
        mmmq = mmmm;                           /* mmmq is to be the mask used to query if a one is under it at this stage, start with old mask */ \
        mmmq &= mmmm^(mmmm<<cccc);             /* divided high part of mask into two equal parts, taking lower part */ \
        tttt = v&mmmq?0:1;                     /* zero if a one under the query mask, one otherwise */ \
        mmmm=mmmq^(tttt*mmmm);                 /* the new mask for next stage is the query mask if a one is under it, otherwise the other half of mmmm */ \
        c+=tttt*cccc;                          /* the right zero counter is incremented by the length of the current interval cccc if a one was not under mask */ \
    }                                          /* note that Anderson's algorithm was incorrect, see also profile comparison in standalone lsb64.c */ \
}
//----------------------------------------------------- list of subroutines -----------------------------------------------------------------------------
// colorgenes1          colour genes specified as input parameters
// colorgenes           colour genes specified at current time point
//.......................................................................................................................................................
// randprob             random event with probability determined by a 32 bit unsigned integer iprob as iprob / 2^32 using RAND128
// selectone_of_2       select one (or none) of two genes based on selection model parameter selection :  returns birth and newgene
// selectone_of_s       select one (or none) of s genes based on selection model parameter selection :  returns birth and newgene
// selectone_nbs        select one of two genes based on pattern of their live 2nd shell neighbours and their genetic encoding
// selectdifft1         select the gene at the single active neighbour position : algorithm could be optimized
// selectdifft2         select the right or left of two genes bunched with least number of empty genes between them
// selectdifft3         select the unique most different (by symmetry) of three live neighbours or first one in canonical rotation
// selectdifft4         select the most central (left) of four live neighbours or first one in canonical rotation
// selectdifft5         select the most central (left) of five live neighbours or first one in canonical rotation
// selectdifft6         select the most central (left) of six live neighbours or first one in canonical rotation
// selectdifft7         select the most central (left) of seven live neighbours or first one in canonical rotation
// selectdifft          select the most central (left) of sum live neighbours or first one in canonical rotation : calls 1 of selectdifft1-7
// disambiguate         disambiguate the cases where the canonical rotation does not uniquely identify a pattern start point : 1 of 8 methods
// testlut              tests LUT calculation via selectdifft: prepares shortcut array lookup
//.......................................................................................................................................................
// patt_hash            hash function for a pattern specified by 4 64-bit (8x8) patterns
// newkey               assign one of free keys kept in a linked list, allocated 1024 at a time (used if hash-key already occupied with different quadtree)
// node_hash            hash function for a node specified by 4 64-bit pointers
// hash_patt_find       find quadtree hash for pattern (leaf of quadtree consists of 4 64bit integers defining a 16x16 pixel array)
// hash_node_find       find quadtree hash for node (node is specified by its four quadrant pointers (64 bit))
// quadimage            construct quadtree for an entire image, reporting if the image has been found previously
//.......................................................................................................................................................
// pack012neighbors     pack all up to 2nd neighbours in single word
// pack0123neighbors    pack all up to 3rd neighbours in single word
// pack49neighbors      fast routine to pack all up to 3rd neighbours in single word : oder of bits dictated by hiearchical assembly
// pack16neighbors      pack 4x4 blocks in single uint64_t word using 16 bits
// pack64neighbors      pack 8x8 blocks in single uint64_t (long) word with 64 bits
// compare_neighbors    compare packed pack neighbours with one given x,y shift of arbitrary size
// compare_all_neighbors compare packed pack neighbours with all nearest neighbour x,y shifts
// packandcompare       pack and compare either all 1-shifted 3-neighbourhoods with t=-1 or chosen (dx,dy,dt) 3-neighbourhoods
//.......................................................................................................................................................
// lab_union            disjoint rank union of equivalence classes returning common root
// label_cell           label a cell (site) in the cellular automata with first pass label of connected component
// label_components     do two-pass fast component labelling with 8-neighbour using Suzuki decision tree, rank union and periodic BCs
// extract_components   extract labelled components to list of subimages embedded in square of side 2^n, each stored in a quadtree hash table
//.......................................................................................................................................................
// update               update gol, golg, golgstats for a single synchronous time step : for selection 0-7 with fixed GoL rule departures in repscheme
// update_lut_sum       update version for gene encoding look up table for totalistic survival and birth (disallowing 0 live neighbour entries) sel 8,9
// update_lut_dist      update version for gene encoding look up table for survival and birth based on corner & edge sums (2*19 states disallowing s=0,1,7,8): sel 10,11
// update_lut_canon_rot update version for gene encoding look up table for canonical rotation survival and birth (2*32 states, disallowing 0,1,7,8 entries) : sel 12,13
// update_lut_2D_sym    update version all different configurations under the standard 2D 4-rotation and 4-reflection symmetries are distinguished: sel 14,15
// update_gol16         update version for 16 parallel gol planes, coupled by a joint gene for all planes : sel 16-19
// update_gol64         update version for 64 parallel gol planes, currently without genetic coupling (routine not yet used) : sel 20-24
// update_gol2match     update version for 2 parallel gol planes, coupled by a gene on one plane by matching of neighborhood : sel 25
//.......................................................................................................................................................
// tracestats           record the current stats in time trace
// get_stats            get the traced statistics from C to python
// countconfigs         count the configs with pthon specified offsets in (x,y,t)
//.......................................................................................................................................................
// get_hist             get the histogram from C to python
// get_activities       get the activity statistics from C to python
//.......................................................................................................................................................
// genelife_update      call update, collect statistics if required and rotate planes
// initialize_planes    initialize periodic sequence of planes to record rolling time window of up to maxPlanes time points (≤8)
// readFile             read file of gol/golg array (32x32) data
// writeFile            write file of gol/golg array (32x32) data
// initialize           initialize simulation parameters and arrays
//.......................................................................................................................................................
// set_colorfunction    set color function integer from GUI for use in patterning and coloring display
// setget_act_ymax      set activity ymax for scaling of activity plot
// set_selectedgene     set selected gene for highlighting from current mouse selection in graphics window
// set_offsets          set offsets for detection of glider structures in display for color function 8
// set_quadrant         set the pair of bits in repscheme (or survivalmask or overwritemask) used for quadrant variation 0-6
// set_randomsoup       change the randomsoup activation for rbackground if nonzero or continual updating of init field with random states and genes : 2,1,0
// set_rbackground      set the backround random live gene input rate per frame and site to rbackground/32768
// set_repscheme_bits   set the two of the repscheme (or survivalmask or overwritemask) bits corresponding to the selected quadrant
// set_repscheme        set repscheme from python
// set_rulemod          set rulemod from python
// set_surviveover      set the two masks for survival and overwrite from python (survivalmask, overwritemask)
// set_vscrolling       set vertical scrolling to track fronts of growth in vertical upwards direction
//.......................................................................................................................................................
// get_log2N            get the current log2N value from C to python
// get_curgol           get current gol array from C to python
// get_curgolg          get current golg array from C to python
// get_acttrace         get current acttrace array C to python
// get_curgolgstats     get current golgstats array C to python
//.......................................................................................................................................................
// cmpfunc              compare gene values as numerical unsigned numbers
// cmpfunc1             compare gene counts in population
// countspecies1        count genes with gene array specified as input parameters
// countspecies         count different genes with genes specified at current time point
// cmpfunc2             compare gene values corresponding to given number index in hash table
// cmpfunc3             compare population counts of hash stored genes
// countspecieshash     count different genes in current population from record of all species that have existed
// activitieshash       calculate array of current activities and update acttrace array of genes in activity plot format
// cmpfunc4             compare birth times of hash stored genes
// cmpfunc5             compare common genealogy level gene values of hash stored genes
// genealogies          calculate and display genealogies
//.......................................................................................................................................................
// get_sorted_popln_act return sorted population and activities (sorted by current population numbers)
// delay                time delay in ms for graphics
// printxy              terminal screen print of array on xterm
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------- colorgenes -------------------------------------------------------------------------------------
void colorgenes1(uint64_t gol[],uint64_t golg[], uint64_t golgstats[], int cgolg[], int NN2) {
    uint64_t gene, gdiff, g2c, mask;
    int ij,k,activity,popcount;
    unsigned int d,d0,d1,d2;
    unsigned int color[3],colormax;
    double rescalecolor;
    static int numones[16]={0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};

    if(colorfunction==0) { // colorfunction based on multiplicative hash
        // see https://stackoverflow.com/questions/6943493/hash-table-with-64-bit-values-as-key/33871291
        for (ij=0; ij<NN2; ij++) {
            if (gol[ij]) {
                gene = golg[ij];
                if (gene == 0ull) gene = 11778L; // random color for gene==0
                // mask = (gene * 11400714819323198549ul) >> (64 - 8);   // hash with optimal prime multiplicator down to 8 bits
                // mask = (gene * 11400714819323198549ul) >> (64 - 32);  // hash with optimal prime multiplicator down to 32 bits
                mask = gene * 11400714819323198549ull;
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
                POPCOUNT64C(gene,d);                                        // assigns number of ones in gene to d
                switch (selection) {
                        case 0 : gdiff=(gene>>40); mask = ((((gdiff&0xff)<<16)+(((gdiff>>8)&0xff)<<8)+(gdiff>>16))<<8) + 0xff; break;  // MSByte red shade, next 8 green, next 8 blue: highest=white
                        case 1 : mask = d==64? 0x0000ffff : ((((63-d)<<18)+(abs((int)d-32)<<10)+(d<<2))<<8) + 0xff; break;  // number of ones determine color gradient blue to red (all 1s red)
                        case 2 :
                        case 3 : d = d & 0x3; mask = d==3 ? 0xf0f0f0ff : ((0xff<<(d<<3))<<8)+0xff; break;  // scissors-stone-well-paper: red-green-blue-white
                        case 4 : mask = d < ncoding ? ((0x3f^d)<<19)+0xff : ((64-d < ncoding) ? ((0x3f^(64-d))<<11)+0xff : 0xf0f0f0ff); break; // near 0 green, near 1 red, others white
                        case 5 : mask = d >= 32 ? ((0x3f^(64-d))<<11)+0xff : ((0x3f^d)<<19)+0xff; break;  //predators green, prey red
                        case 6 : g2c = (1ull<<ncoding)-1ull;gdiff = gene^g2c; POPCOUNT64C(gdiff,d2);
                                 mask = d<d2 ? (d<<26)+0xff : (d2<<10)+0xff; break;
                        case 7 : g2c = (gene>>8)&((1ull<<ncoding)-1ull);
                                 gdiff = gene&0xff;POPCOUNT64C(gdiff,d2);d = d2>7? 7 : d2;
                                 mask = g2c ? 0xf0f0f0ff : ((0x1f+(d<<5))<<8)+(((gdiff>>4)&0xf)<<27)+(((gdiff&0xf)<<4)<<16)+0xff; break;
                        case 8 :                                            // colors according to nr of non zero LUT entries d from blue to red in increasing d
                                 if (ncoding==1) {gdiff=gene&0xffff;POPCOUNT64C(gdiff,d);mask = d==16 ? 0x0000ffff : ((((15-d)<<20)+(abs((int)d-8)<<12)+(d<<4))<<8) + 0xff;}
                                 else if (ncoding==2) {PATTERN2_32(gene,0x3,d);mask = d==16 ? 0x00ffffff : ((((15-d)<<20)+(abs((int)d-8)<<12)+(d<<4))<<8) + 0xff;}
                                 else {PATTERN4(gene,0xf,d);mask = d==16 ? 0x00ffffff : ((((15-d)<<20)+(abs((int)d-8)<<12)+(d<<4))<<8) + 0xff;} break;
                        case 9 : PATTERN4(gene,0x0,d0);PATTERN4(~gene&0x8888888888888888ull,0x8,d1);PATTERN4(gene&0x8888888888888888ull,0x8,d2);
                                 mask = (d2<<26)+((d1-d0)<<10)+0xff; break;
                        case 10: POPCOUNT64C(gene&0x7ffffull,d1); POPCOUNT64C((gene>>32)&0x7ffffull,d2);
                                 mask = (d2<<27)+(d1<<11)+0xff; break;
                        case 11: PATTERN8(gene,0x27,d0);d0=d0>1?2:d0;PATTERN8(gene,0x3f,d1);d1=d1>1?2:d1;d0=d0+d1;d0=d0>2?3:d0; // saturated copy nrs for survival GoL LUTs
                                 PATTERN8(gene,0xbf,d1);d1=d1>2?3:d1; PATTERN8(gene&0x8080808080808080ull,0x80,d2);d2=d2>7?7:d2; // saturated copy nrs for birth GoL and birth LUTs
                                 mask =(((d0<<22)+(d1<<14)+(d2<<5)) << 8) + 0xff; break;
                        case 12: POPCOUNT64C(gene&0xffffffffull,d1); POPCOUNT64C((gene>>32)&0xffffffffull,d2);
                                 mask = (d2<<27)+(d1<<11)+0xff; break;
                        case 13: PATTERN8(gene,0x2f,d0);d0=d0>1?2:d0;PATTERN8(gene,0x3f,d1);d1=d1>1?2:d1;d0=d0+d1;d0=d0>2?3:d0; // saturated copy nrs for survival GoL LUTs
                                 PATTERN8(gene,0xbf,d1);d1=d1>2?3:d1; PATTERN8(gene&0x8080808080808080ull,0x80,d2);d2=d2>7?7:d2; // saturated copy nrs for birth GoL and birth LUTs
                                 mask =(((d0<<22)+(d1<<14)+(d2<<5)) << 8) + 0xff; break;  //Check NYI
                        case 14: POPCOUNT64C(gene&0xffffffffull,d1); POPCOUNT64C((gene>>32)&0xffffffffull,d2);
                                 mask = (d2<<27)+(d1<<11)+0xff; break;
                        case 15: mask = 0xffffffff;  break; // NYI
                        case 16:
                        case 17: for (d=0,mask=0;d<16;d++) {d2=(gene>>(d<<2))&0xf;d2=(d2!=d)?1:0;mask+=(d2*0x2a)<<((d%3)<<3);}
                                 mask = (mask<<8)+0x7f7f7fff;break;  // 0x2a is 42, approx 1/6 (ie <16/3) of 255
                        case 18:
                        case 19: for (d=0,mask=0;d<16;d++) {d2=(gene>>(d<<2))&0xf;d2=numones[d2];mask+=(d2*0xb)<<((d%3)<<3);}
                                 mask = (mask<<8)+0x7f7f7fff;break;  // 0x2a is 42, approx 1/6 (ie <16/3) of 255
                        case 20:
                        case 21:
                        case 22:
                        case 23: d2= (d>>5) ? 5*32 + d-32 : ((d>>4) ? 4*32 + (d-16)*2
                                                         : ((d>>3) ? 3*32 + (d-8)*4
                                                         : ((d>>2) ? 2*32 + (d-4)*8
                                                         : ((d>>1) ? 32 + (d-2)*16
                                                         : d*32))));
                                mask = gene * 11400714819323198549ul;mask = mask >> (64 - 16);
                                mask = ((d2+62)<<8)+(mask<<16)+0xff;break;
                        case 24:
                        case 25: mask = ((gene>>40)<<8)+0xff; break;
                        default  : mask = ((d+(d<<6)+(d<<12)+(d<<18))<<8) + 0xff;
                }
                if(colorfunction==2) {
                    if(golgstats[ij]&F_nongolchg) mask = 0x00ffffff;        // color states changed by non GoL rule yellow
                    if(selection>=20) {                               // color as superposition of multiplane gol states
                        if (displayoneplane==64) {
                            POPCOUNT64C(gol[ij],d);
                            d2= (d>>5) ? 5*32 + d-32 : ((d>>4) ? 4*32 + (d-16)*2
                                                     : ((d>>3) ? 3*32 + (d-8)*4
                                                     : ((d>>2) ? 2*32 + (d-4)*8
                                                     : ((d>>1) ? 32 + (d-2)*16
                                                     : d*32))));
                            mask = gol[ij] * 11400714819323198549ul;mask = mask >> (64 - 16); // 16 random bits based on gol state
                            mask = ((d2+62)<<24)+(mask<<16)+0xff;
                        }
                        else mask = ((unsigned int)((gol[ij]>>displayoneplane)&1ull))*((displayoneplane<<26)+((63-displayoneplane)<<18)+0xff);
                    }
                    else if(selection>=16) {                               // color as superposition of multiplane gol states
                        POPCOUNT64C(gol[ij]&displayplanes,d);
                        d2= (d>>3) ? 3*32 + (d-8)*4
                                    : ((d>>2) ? 2*32 + (d-4)*8
                                    : ((d>>1) ? 32 + (d-2)*16
                                    : d*32));
                        mask = (gol[ij]&displayplanes) * 11400714819323198549ul;mask = mask >> (64 - 16);
                        mask = ((d2+62)<<8)+(mask<<16)+0xff;
                    }
                    else if(selection>=14) {                                     // color as superposition of biplane gol states
                        if (gol[ij]>>1) {
                            gene = golg[ij];
                            if (gene == 0ull) gene = 11778L; // random color for gene==0
                            mask = gene * 11400714819323198549ul;
                            mask = mask >> (64 - 8);   // hash with optimal prime multiplicator down to 8 bits
                        }
                        else mask = 0;
                        mask = ((0xff*(gol[ij]&1ull))<<24)+((0xff*((gol[ij]>>1)&1ull))<<16)+(mask<<8)+0xff;
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
    else if(colorfunction==4){                                  //activities
        int popmax = 100;                                      // need to bring this parameter up to python
        for (ij=0; ij<NN2; ij++) {
            gene=acttrace[ij];
            if (gene == rootgene) mask = 0x3f3f3fff;            // grey color for background, all root genes
            else {
                if (gene == 0ull) gene = 11778L;                // random color for gene==0
                mask = gene * 11400714819323198549ul;
                mask = mask >> (64 - 32);                       // hash with optimal prime multiplicator down to 32 bits
                mask |= 0x808080ff; // ensure brighter color at risk of improbable redundancy, make alpha opaque
                if(popmax) {
                    popcount=0;
                    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) popcount = genedataptr->popcount;
                    else fprintf(stderr,"gene not found in colorfunction for activities\n");
                    if(popcount>popmax) popcount=popmax;
                    colormax=0;
                    for(d=0;d<3;d++) if((color[d]=( (mask>>(8+(d<<3))) & 0xff))>colormax) colormax=color[d];
                    rescalecolor=(log((double)popcount)/log((double)popmax))*((double)0xff/(double)colormax);
                    for(d=0;d<3;d++) color[d]=(unsigned int) (((double) color[d])*rescalecolor);
                    for(d=0,mask=0xff;d<3;d++) mask |= color[d]<<((d<<3)+8);
                }
            }
            cgolg[ij]= (int) mask;
        }
    }
    else if(colorfunction<8){                                       //genealogies
        for (ij=0; ij<NN2; ij++) {
            gene=genealogytrace[ij];
            activity = 0;
            if (gene == selectedgene) mask = 0xffffffff; else
            if (gene == rootgene) mask = 0x000000ff;                // black color for root
            else {
                if(colorfunction==7) {
                    activity = 0;
                    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) activity = genedataptr->activity;
                    else fprintf(stderr,"gene not found in colorfunction for genealogy\n");
                }
                if (gene == 0ull) gene = 11778L;                    // random color for gene==0
                mask = gene * 11400714819323198549ul;
                mask = mask >> (64 - 32);                           // hash with optimal prime multiplicator down to 32 bits
                mask |= 0x808080ff;                                 // ensure brighter color at risk of improbable redundancy, make alpha opaque
                if(colorfunction==7) {                              // rescale color brightness by activity/activitymax
                    colormax=0;
                    for(d=0;d<3;d++) if((color[d]=( (mask>>(8+(d<<3))) & 0xff))>colormax) colormax=color[d];
                    rescalecolor=((double)(activity*0xff))/((double)(activitymax*colormax));         // integer version doesn't work
                    for(d=0;d<3;d++) color[d]=(unsigned int) (((double) color[d])*rescalecolor);     // rescale colors by activity/activitymax
                    for(d=0,mask=0xff;d<3;d++) mask |= color[d]<<((d<<3)+8);
                }
            }
            cgolg[ij]=(int) mask;
        }
    }
    else if(colorfunction==8) {         // colorfunction based on packed bit pattern with multiplicative hash, compare with offset
        for (ij=0; ij<NN2; ij++) {
                gene = golmix[ij];
                if(offdx==0 && offdy==0 && offdt==0) {
                    for (mask=0,k=0;k<8;k++) {
                        d1 = (gene>>(k<<3))&0xff;
                        d1 = (d1 > 63) ? 0 : (63-d1);
                        d1 = (d1 < 48) ? 0 : d1-48;                 // 0 to 15 : perfect match is 15  (4 bits)
                        d1 = (d1==0xf) ? 0x1f : d1;                 // perfect match separated to value 31 (5 bits) for better contrast
                        if(k<3) mask+=d1<<(3+(k<<3));               // perfect match has full intensity colour
                        else if (k==3 && d1==0x1f) mask = (d1<<3)+(d1<<11)+(d1<<19); // the fourth channel has white colour : no others shown
                        else if (k<7 && d1==0x1f) mask+= (d1<<(3+((k-4)<<3)))+(d1<<(3+((k<6?k-3:0)<<3))); // mixed colours for NE SE SW
                        else if(d1==0x1f) mask+= (d1<<3)+(d1<<10)+(d1<<18); // mixed colour for NW
                    }
                    mask = (mask<<8)+0xff;
                }
                else {
                   POPCOUNT64C(gene,d);                    // assigns number of ones in gene to d. These 3 lines version for one offset comparison
                    d=(d==64)?0:63-d;
                    mask = (d==63) ? 0xffffffff : ((((d&3)<<22)+(((d>>2)&3)<<14)+(((d>>4)&3)<<6))<<8) + 0xff;
                }

                cgolg[ij] = (int) mask;
        }
    }
    else if(colorfunction==9) {         // colorfunction based on unique labelling of separate components in image
        for (ij=0; ij<NN2; ij++) {
            if (gol[ij]) {
                gene = (uint64_t) label[ij];
                mask = gene * 11400714819323198549ull;
                mask = mask >> (64 - 32);   // hash with optimal prime multiplicator down to 32 bits
                mask |= 0x808080ff; // ensure brighter color at risk of improbable redundancy, make alpha opaque
                cgolg[ij] = (int) mask;
            }
            else cgolg[ij] = 0;
        }
    }
}
//.......................................................................................................................................................
void colorgenes(int cgolg[], int NN2) {
    colorgenes1(gol, golg, golgstats, cgolg, NN2);
}
//-------------------------------------------------------- randprob ------------------------------------------------------------
extern inline uint64_t randprob(unsigned int uprob, unsigned int randnr) {
    return(randnr < uprob ? 1ull : 0ull);
}
//------------------------------------------------------- selectone ------------------------------------------------------------
extern inline void selectone_of_2(int s, uint64_t nb2i, int nb[], uint64_t golg[], uint64_t * birth, uint64_t *newgene, unsigned int kch) {
// birth is returned 1 if ancestors satisfy selection condition. Selection of which of two genes to copy is newgene. Non-random result.
    unsigned int k,d0,d1,d2,d3,dd,swap;                  // number of ones in various gene combinations
    uint64_t livegenes[2],gdiff,gdiff0,gdiff1; // various gene combinations
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
            *birth = (d0^d1) ? 1ull: 0ull;               // birth condition is two genes different in number of ones
            *newgene= (d0>d1) ? livegenes[0] : livegenes[1];
            break;
        case 2:                                          // scissors-stone-well-paper game on number ones mod 4
                                                         // scissors 0 stone 1 well 2 paper 3
                                                         // exception to numerical order: sc>pa
            d0=d0&0x3;d1=d1&0x3;
            *birth = (d0^d1) ? 1ull: 0ull;               // birth if 2 genes differ mod 4 in number of ones
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
            *birth = ((d0^d1)==1ull) ? 1ull: 0ull;       // birth if 2 genes differ by 1 mod 4 in number of ones
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
        case 4:                                          // birth if 2 genes cooperate : closer to all 0, all 1 targets than ncoding and closer to each other than 64-ncoding
            gdiff=livegenes[0]^livegenes[1];
            POPCOUNT64C(gdiff,dd);
            *birth = dd>0 && dd<(64-ncoding) && ((d0<ncoding && d1>64-ncoding) || (d1<ncoding && d0>64-ncoding))  ? 1ull: 0ull; // birth if 2 genes close enough to targets
            if (d0<ncoding) {if(d0>64-d1) swap=1;else swap=0;}                  // need 64-d1 != d0 to avoid asymmetry in direction               && (64-d1 != d0)
            else {if(64-d0>=d1) swap=1; else swap=0;}
            *newgene= livegenes[swap];
            break;
        case 5:                                          // predator prey model: prey-prey evolves to all 0, predator to complement of prey
            gdiff=livegenes[0]^livegenes[1];
            gdiff1=livegenes[0]^(~livegenes[1]);
            POPCOUNT64C(gdiff1,dd);
            prey = (d0<32) || (d1<32);                   // prey present : newgene is one with less ones, 1 prey : predator wins
            prey2 = (d0<32) && (d1<32);                  // 2 prey : newgene is one with less ones, 1 prey : predator wins
            // *birth = (gdiff && prey && dd<ncoding) ? 1ull: 0ull;  // birth if different and >=1 prey and close enough match)
            *birth = (gdiff && prey2) || (prey && (!prey2) && (dd<ncoding)) ? 1ull: 0ull; // birth if different and >=1 prey and close enough match)
            *newgene= (prey2 ? ((d0<d1) ? livegenes[0] : livegenes[1]) : ((d0<32) ? livegenes[1] : livegenes[0]));
            break;
        case 6:                                         // birth if 2 genes differently functional (Not Yet Working)
                                                                            // 1st target is WLOG the all 0 sequence, d0 and d1 gives distances of two genes
            if(ncoding<64) gene2centre = (1ull<<ncoding)-1ull;              // first ncoding 1s in this sequence
            else gene2centre = ~0ull;                                       // otherwise 2nd target is all 64 ones
            gdiff  = livegenes[0]^livegenes[1];                             // difference between two genes
            gdiff0 = livegenes[0]^gene2centre;                              // difference of gene 0 to 2nd target
            gdiff1 = livegenes[1]^gene2centre;                              // difference of gene 1 to 2nd target
            POPCOUNT64C(gdiff,dd);                                          // dd is distance between genes
            POPCOUNT64C(gdiff0,d2);                                         // d2 is distance of gene 0 from 2nd target
            POPCOUNT64C(gdiff1,d3);                                         // d3 is distance of gene 1 from 2nd target
            g0011 = d0<dd && d3<dd;                      // gene 0 closer to 1st target and gene 1 closer to 2nd target than they are from each other
            g0110 = d2<dd && d1<dd;                      // gene 0 closer to 2nd target and gene 1 closer to 1st target than they are from each other
            *birth = (g0011 && (d0!=d3)) != (g0110 && (d2!=d1))  ? 1ull: 0ull; // birth if 2 genes closer to two different targets than to each other
            *newgene= (g0011 && (d0!=d3)) ? ((d0<d3) ? livegenes[0] : livegenes[1]) : ((d2<d1) ? livegenes[0] : livegenes[1]);
            break;
        case 7:                                          // neutral selection but selective birth only occurs if two chosen sequences are same (different) (NB uses RNG)
                                                         // note that birth is already suppressed for 3 identical live nbs unless enforcebirth-3 bit on
                                                         // we want to use this routine only for two live nbs and only when genes not same
            *birth = livegenes[0]^livegenes[1] ? 1ull: 0ull;
            *newgene = golg[nb[kch]];
            break;
        // cases 8+: this subroutine is not used
        default:
            fprintf(stderr,"Error: two live gene fitness value %d is not implemented\n",selection);
            exit(1);
    }
}
//.......................................................................................................................................................
extern inline int selectone_of_s(int s, uint64_t nb1i, int nb[], uint64_t golg[], uint64_t * birth, uint64_t *newgene, uint64_t *nbmask) {
// result is number of equally fit best neighbours that could be the ancestor (0 if no birth allowed, 1 if unique) and list of these neighbour indices
// birth is returned 1 if ancestors satisfy selection condition. Selection of which of genes to copy is newgene. Non-random result.
    unsigned int k,nbest;                                 // index for neighbours and number in best fitness category
    unsigned int d[8],dS,dB,d0,d1,d2;                     // number of ones in various gene combinations
    uint64_t livegenes[8],gdiff,extrval,bestnbmask;
    unsigned int repselect = (repscheme>>4)&0x7;

    if(s==0) {*birth = 1ull;  *newgene = genegol[selection-8];return(1);}

    for(k=0;k<s;k++) {
        livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
        POPCOUNT64C(livegenes[k],d[k]);
    }

    switch (repselect) {
        case 0:                                          // integer value of sequence as fitness : smallest or largest
        case 1:
            if(repselect&0x1) for(extrval= 0ull,k=0;k<s;k++) extrval = (livegenes[k] >= extrval ? livegenes[k] : extrval); // find value of fittest gene max value
            else              for(extrval=~0ull,k=0;k<s;k++) extrval = (livegenes[k] <= extrval ? livegenes[k] : extrval); // find value of fittest gene min value
            for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= ((livegenes[k]==extrval) ? 1ull<< (k+0*nbest++) : 0ull); // find set of genes with equal best value
            *birth = (nbest>0) ? 1ull: 0ull;             // birth condition may include later that genes not all same
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;   // execute loop until first optimal live neighbour found at k (a one)
            *newgene = livegenes[k&0x7];                 // choose first of selected set to replicate (can make positional dependent choice instead externally)
            break;
        case 2:
        case 3:                                          // number of ones in sequence as fitness
            if(repselect&0x1) for(extrval= 0ull,k=0;k<s;k++) extrval = (d[k]>= extrval ? d[k] : extrval); // find value of fittest gene max nr ones
            else              for(extrval=~0ull,k=0;k<s;k++) extrval = (d[k]<= extrval ? d[k] : extrval); // find value of fittest gene min nr ones
            for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= (d[k]==extrval ? 1ull<<(k+0*nbest++) : 0ull); // find set of genes with equal best value
            *birth = ((nbest>0) ? 1ull: 0ull);           // birth condition may include later that genes not all same
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            *newgene = livegenes[k&0x7];                 // choose first of selected set to replicate (can make positional dependent choice instead externally)
            break;
        case 4:                                          // neutral selection
            for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= 1ull<<(k+0*nbest++); // find set of genes with equal best value : in this case all of them
            *birth = 1ull;                               // birth condition is always true
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            *newgene = livegenes[k&0x7];                 // choose first of selected set to replicate (can make positional dependent choice instead externally)
            break;
        case 5:                                          // neutral selection but birth only occurs if some sequences are different
            for(gdiff=0ull,k=1;k<s;k++) if ((gdiff=livegenes[0]^livegenes[k])) break; // test whether all genes the same
            if (gdiff) for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= 1ull<<(k+0*nbest++); // find set of genes with equal best value : in this case all of them
            else { nbest = 0; bestnbmask = 0ull;}
            *birth = gdiff ? 1ull : 0ull;                // birth condition is genes not all same
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            if (k==s) {k=0;*birth = 0ull;}               // avoid analyze error, and be on safe side
            *newgene = livegenes[k&0x7];                 // choose first of selected set to replicate (can make positional dependent choice instead externally)
            break;
        case 6:
        case 7:
          switch(selection) {
            case 8:                                          // totalistic lut penalty of gene in fixed length encoding : first 8 bits survival, next 8 bits birth
                for(k=0;k<s;k++) {
                    livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    POPCOUNT64C(livegenes[k]&0xffull,dS);
                    POPCOUNT64C((livegenes[k]>>8)&0xffull,dB);
                    d[k]=24-dS-(dB<<1);
                }
                break;
            case 9:                                          // totalistic lut penalty of gene in variable length encoding : 4 bit patterns S 0xxx  B 1xxx
                for(k=0;k<s;k++) {
                    livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    PATTERN4(livegenes[k],0x0,d0)
                    PATTERN4(~livegenes[k]&0x8888888888888888ull,0x8,d1)
                    PATTERN4(livegenes[k]&0x8888888888888888ull,0x8,d2)
                    d[k]=48-(d1-d0)-2*d2;
                }
                break;
            case 10:
                for(k=0;k<s;k++) {
                    livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    POPCOUNT64C(livegenes[k]&0x7fffffull,dS);        // 23 coding bits for survival
                    POPCOUNT64C((livegenes[k]>>32)&0x7fffffull,dB);  // 23 coding bits for birth
                    d[k]=57-dS-(dB<<1);
                }
                break;
            case 11:                                        // dist. lut penalty of gene in variable length encoding : 8 bit patterns S 0xxxrrrr  B 1xxxrrrr
                for(k=0;k<s;k++) {
                    livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    PATTERN8(livegenes[k],0x00,d0)
                    PATTERN8(~livegenes[k]&0x8080808080808080ull,0x80,d1)
                    PATTERN8(livegenes[k]&0x8080808080808080ull,0x80,d2)
                    d[k]=24-(d1-d0)-2*d2;
                }
                break;
            case 12:                                        // canon. lut penalty in fixed length encoding 32 bit survival 32 bit birth
            case 14:
                 for(k=0;k<s;k++) {
                    livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    POPCOUNT64C(livegenes[k]&0xffffffffull,dS);        // 32 coding bits for survival
                    POPCOUNT64C((livegenes[k]>>32)&0xffffffffull,dB);  // 32 coding bits for birth
                    d[k]=96-dS-(dB<<1);
                }
                break;
            case 13:                                        // canon. lut penalty of gene in variable length encoding : 8 bit patterns S 0xxxwrrr  B 1xxxwrrr
                for(k=0;k<s;k++) {                          // could be made more accurate by accounting for upper pair bit states for rest of 5-subset too
                    livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    PATTERN8(livegenes[k]&0xffffffffffff,0x00,d0)
                    PATTERN8(~livegenes[k]&0x808080808080ull,0x80,d1)
                    PATTERN8(livegenes[k]&0x808080808080ull,0x80,d2)
                    d[k]=24-(d1-d0)-2*d2;
                }
                break;
            case 15:                                        // dist. or canon. lut penalty of gene in variable length encoding : 8 bit patterns S 0xxxwwrr  B 1xxxwwrr
                for(k=0;k<s;k++) {                          // could be made more accurate by accounting for upper quartet bit states for rest of 6-subset too
                    livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    PATTERN8(livegenes[k]&0xffffffffff,0x00,d0)
                    PATTERN8(~livegenes[k]&0x8080808080ull,0x80,d1)
                    PATTERN8(livegenes[k]&0x8080808080ull,0x80,d2)
                    d[k]=24-(d1-d0)-2*d2;
                }
                break;
            default:
                fprintf(stderr,"Error: selectone_of_s for selection %d is not implemented\n",selection);
                exit(1);
          }
          for(extrval=0ull,k=0;k<s;k++) extrval = d[k]>= extrval ? d[k] : extrval; // find value of fittest gene
          for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= (d[k]==extrval) ? 1ull<<(k+0*nbest++) : 0ull; // find set of genes with equal best value
          *birth = 1ull;                               // birth condition is always met since extrval set is always > 0
          for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
          *newgene = livegenes[k&0x7];                 // choose first of selected set to replicate (can make positional dependent choice instead externally)
          break;
        default:
            fprintf(stderr,"Error: s = %d live gene repselect %d is not implemented\n",s,repselect);
            exit(1);
    }
    for (*nbmask=0ull,k=0;k<s;k++)
        *nbmask |= ((bestnbmask>>k)&0x1ull)<<((nb1i>>(k<<2))&0x7);
    return(nbest);
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
//------------------------------------------------------- selectdifft0 -------------------------------------------------------------------------------------
extern inline unsigned int selectdifft0(uint64_t nbmask, int *crot, int *kodd) {
    *kodd = 0;
    *crot = 0;
    return(0);                                               // replication of live nb in bit 0 of canonical rotation
}
//------------------------------------------------------- selectdifft1 -------------------------------------------------------------------------------------
extern inline unsigned int selectdifft1(uint64_t nbmask, int *crot, int *kodd) {
// selection based on canonical rotation
    int k,kmin;
    uint64_t nbmaskr, nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {       // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);    // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                            // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                              // neighbor mask rotate min is current rotation
            kmin = k;                                        // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;
    *crot = 0;
    return(kmin);                                            // replication of live nb in bit 0 of canonical rotation
}
//.......................................................................................................................................................
extern inline unsigned int selectdifft2(uint64_t nbmask, int *crot, int *kodd) {
// selection based on canonical rotation to bunched pair, choose clockwise or anti-clockwise one (for R_2_canonical_nb)
    int k,kmin;
    uint64_t nbmaskr, nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {       // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);    // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                            // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                              // neighbor mask rotate min is current rotation
            kmin = k;                                        // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;
    switch (nbmaskrm) {                                        //              x03    x05    x09    x11
        case 0x03ull : k = 1; *crot = 0; break;                // 00000011    |01.|  <-
        case 0x05ull : k = 2; *crot = 1; break;                // 00000101    |...|  |0.2|  <-
        case 0x09ull : k = 3; *crot = 2; break;                // 00001001    |...|  |...|  |0..|   <-
        case 0x11ull : k = 4; *crot = 3; break;                // 00010001           |...|  |..3|  |0..|   <-
        default  : {                                           //                           |...|  |...|
                                                               //                                  |..4|
            fprintf(stderr,"Error in canonical rotation for two live neighbours nbmaskrm = %llx for mask %llx\n",nbmaskrm,nbmask); k = 0;
        } //default case
    } //switch

    if (canonical) return(kmin);                                      // replication of live neigbour in bit 0 of canonical rotation
    else return((kmin+k)&0x7);                                        // rotate unique nb k left (kmin) back to orig nb pat
}
//.......................................................................................................................................................
extern inline unsigned int selectdifft3(uint64_t nbmask, int *crot, int *kodd) {
    unsigned int k,kmin;
    uint64_t nbmaskr,nbmaskrm;
    
    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {                // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 1ull)<<7) | (nbmaskr>>1);               // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                                     // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                       // neighbor mask rotate min is current rotation
            kmin = k;                                                 // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;                                                                     // replication of live neighbour in most different position
    switch (nbmaskrm) {                              //              x07    x0b    x0d    x13    x15    x19    x25
        case 0x07ull : k = 1; *crot = 0; break;      // 00000111    |012|  <-
        case 0x0bull : k = 0; *crot = 1; break;      // 00001011    |...|  |01.|  <-
        case 0x0dull : k = 3; *crot = 2; break;      // 00001101    |...|  |..3|  |0.2|   <-
        case 0x13ull : k = 1; *crot = 3; break;      // 00010011           |...|  |..3|  |01.|   <-
        case 0x15ull : k = 2; *crot = 4; break;      // 00010101                  |...|  |...|  |0.2|   <-
        case 0x19ull : k = 0; *crot = 5; break;      // 00011001                         |..4|  |...|  |0..|   <-
        case 0x25ull : k = 5; *crot = 6; break;      // 00100101                                |..4|  |..3|  |0.2|  <-
        default  : {                                 //                                                |..4|  |...|
                                                     //                                                       |.5.|
            fprintf(stderr,"Error in canonical rotation for three live neighbours nbmaskrm = %llx for mask %llx\n",nbmaskrm,nbmask); k = 0;
        } //default case
    } //switch
    if (canonical) return(kmin);                                      // replication of live neigbour in bit 0 of canonical rotation
    else return((kmin+k)&0x7);                                        // rotate unique nb k left (kmin) back to orig nb pat

}
//.......................................................................................................................................................
extern inline unsigned int selectdifft4(uint64_t nbmask, int *crot, int *kodd) {
    int k,kmin;
    uint64_t nbmaskr,nbmaskrm;
    
    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {                // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);             // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                                     // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                       // neighbor mask rotate min is current rotation
            kmin = k;                                                 // no of times rotated to right
        }
    }
    *kodd = kmin&0x1;                                                                 // replication of live neighbour in most central position (left disambiguation)
    switch (nbmaskrm) {                             //              x07f    x17    x1b    x1d    x27    x2b    x2d   x33   x35    x55
        case 0x0full : k = 1; *crot = 0; break;     // 00001111    |012|  <-
        case 0x17ull : k = 2; *crot = 1; break;     // 00010111    |..3|  |012|  <-
        case 0x1bull : k = 1; *crot = 2; break;     // 00011011    |...|  |...|  |01.|   <-
        case 0x1dull : k = 2; *crot = 3; break;     // 00011101           |..4|  |..3|  |0.2|   <-
        case 0x27ull : k = 2; *crot = 4; break;     // 00100111                  |..4|  |..3|  |012|   <-
        case 0x2bull : k = 3; *crot = 5; break;     // 00101011                         |..4|  |...|  |01.|   <-
        case 0x2dull : k = 2; *crot = 6; break;     // 00101101                                |.5.|  |..3|  |0.2|  <-
        case 0x33ull : k = 1; *crot = 7; break;     // 00110011                                       |.5.|  |..3|  |01.|  <-
        case 0x35ull : k = 2; *crot = 8; break;     // 00110101                                              |.5.|  |...|  |0.2|  <-
        case 0x55ull : k = 2; *crot = 9; break;     // 01010101                                                     |.54|  |...|  |0.2|  <-
        default  : {                                //                                                                     |.54|  |...|
                                                    //                                                                            |6.4|
            fprintf(stderr,"Error in canonical rotation for four live neighbours nbmaskrm = %llx for mask %llx\n",nbmaskrm,nbmask); k = 0;
        } //default case
    } //switch
    if (canonical) return(kmin);                                      // replication of live neigbour in bit 0 of canonical rotation
    else return((kmin+k)&0x7);                                        // rotate unique nb k left (kmin) back to orig nb pat
}
//.......................................................................................................................................................
extern inline unsigned int selectdifft5(uint64_t nbmask, int *crot, int *kodd) {
    unsigned int k,kmin;
    uint64_t nbmaskr,nbmaskrm;
    
    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {                // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 1ull)<<7) | (nbmaskr>>1);               // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                                     // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                       // neighbor mask rotate min is current rotation
            kmin = k;                                                 // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;
    switch (nbmaskrm) {                              //              x1f    x2f    x3d    x3b    x57    x37    x5b
        case 0x1full : k = 2; *crot = 0; break;      // 00011111    |012|  <-
        case 0x2full : k = 2; *crot = 1; break;      // 00101111    |..3|  |012|  <-
        case 0x3dull : k = 3; *crot = 2; break;      // 00111101    |..4|  |..3|  |0.2|   <-
        case 0x3bull : k = 3; *crot = 3; break;      // 00111011           |.5.|  |..3|  |01.|   <-
        case 0x57ull : k = 2; *crot = 4; break;      // 01010111                  |.54|  |..3|  |012|   <-
        case 0x37ull : k = 2; *crot = 5; break;      // 00110111                         |.54|  |...|  |012|   <-
        case 0x5bull : k = 3; *crot = 6; break;      // 01011011                                |6.4|  |...|  |01.|  <-
        default  : {                                 //                                                |.54|  |..3|
                                                     //                                                       |6.4|
            fprintf(stderr,"Error in canonical rotation for five live neighbours nbmaskrm = %llx for mask %llx\n",nbmaskrm,nbmask); k = 0;
        } //default case
    } //switch
    if (canonical) return(kmin);                                      // replication of live neigbour in bit 0 of canonical rotation
    else return((kmin+k)&0x7);                                        // rotate unique nb k left (kmin) back to orig nb pat
}
//.......................................................................................................................................................
extern inline unsigned int selectdifft6(uint64_t nbmask, int *crot, int *kodd) {
// selection based on canonical rotation to bunched pair, choose clockwise or anti-clockwise one (for R_2_canonical_nb)
    int k,kmin;
    uint64_t nbmaskr, nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {       // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);    // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                            // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                              // neighbor mask rotate min is current rotation
            kmin = k;                                        // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;
    switch (nbmaskrm) {                                        //              x3f    x5f    x6f    x77
        case 0x3full : k = 2; *crot = 0; break;                // 00111111    |012|  <-
        case 0x5full : k = 3; *crot = 1; break;                // 01011111    |..3|  |012|  <-
        case 0x6full : k = 3; *crot = 2; break;                // 01101111    |.54|  |..3|  |012|   <-
        case 0x77ull : k = 4; *crot = 3; break;                // 01110111           |6.4|  |..3|  |012|   <-
        default  : {                                           //                           |65.|  |...|
                                                                   //                                  |654|
            fprintf(stderr,"Error in canonical rotation for six live neighbours nbmaskrm = %llx for mask %llx\n",nbmaskrm,nbmask); k = 0;
        } //default case
    } //switch
    if (canonical) return(kmin);                                      // replication of live neigbour in bit 0 of canonical rotation
    else return((kmin+k)&0x7);                                        // rotate unique nb k left (kmin) back to orig nb pat
}
//...................................................................................................................................................
extern inline unsigned int selectdifft7(uint64_t nbmask, int *crot, int *kodd) {
// selection based on canonical rotation
    int k,kmin;
    uint64_t nbmaskr, nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {       // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);    // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                            // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                              // neighbor mask rotate min is current rotation
            kmin = k;                                        // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;
    *crot = 0;
    if (canonical) return(kmin);                            // replication of live nb in bit 0 of canonical rotation
    else  return((kmin+3)&0x7);                             // replication of live nb in bit 3 (middle) of canonical rotation
}
//...................................................................................................................................................
extern inline unsigned int selectdifft(int sum, uint64_t nbmask, int *crot, int *kodd, int *nsame) {
        int kch;
        *nsame = 0;
    
        switch(sum) {
                    case 0:  return(selectdifft0(nbmask, crot, kodd));
                    case 1:  return(selectdifft1(nbmask, crot, kodd));
                    case 2:  kch=selectdifft2(nbmask, crot, kodd);
                             if (*crot==3) *nsame = 2;
                             return(kch);
                    case 3:  return(selectdifft3(nbmask, crot, kodd));
                    case 4:  kch=selectdifft4(nbmask, crot, kodd);
                             if (*crot==2) *nsame = 2;
                             else if (*crot==9) *nsame = 4;
                             return(kch);
                    case 5:  return(selectdifft5(nbmask, crot, kodd));
                    case 6:  kch=selectdifft6(nbmask, crot, kodd);
                             if (*crot==3) *nsame = 2;
                             return(kch);
                    case 7:  return(selectdifft7(nbmask, crot, kodd));
                    default: return(0);
        }
}
//...................................................................................................................................................
extern inline uint64_t disambiguate(unsigned int kch, uint64_t nb1i, int nb[], uint64_t golg[], int nsame, uint64_t *birth, uint64_t randnr) {
    uint64_t gene,newgene;
    unsigned int d,dmin;
    int k;
    
    switch ((repscheme>>8)&0x7) {
        case 0:  kch += ((nsame-1)&(randnr>>32))<< (nsame ==4 ? 1 : 2);                  // random choice
                 kch &=0x7; return( golg[nb[(nb1i>>(kch<<2))&0x7]]);
        case 1:  return( golg[nb[(nb1i>>(kch<<2))&0x7]]);                                // ignore asymmetry issue, continue regardless;
        case 2:  birth = 0; return(0ull);                                                // abandom birth attempt
        case 3:  for (newgene=~0ull,k=0;k<nsame;k++) {                                   // choose minimum value gene
                     kch+=k*(nsame==4 ? 2 : 4);
                     gene=golg[nb[(nb1i>>(kch<<2))&0x7]];
                     if (gene<newgene) newgene = gene;
                 }; return(newgene);
        case 4:  for (newgene=~0ull,dmin=64+1,k=0;k<nsame;k++) {                         // choose least nr ones gene or minimum value if same
                     kch+=k*(nsame==4 ? 2 : 4);
                     gene=golg[nb[(nb1i>>(kch<<2))&0x7]];
                     POPCOUNT64C(gene,d);
                     if(d<dmin) newgene = gene;
                     else if (d==dmin && newgene < gene) newgene = gene;
                 }; return(newgene);
        case 5:  for (newgene=~0ull,k=0;k<nsame;k++) {                                   // choose AND of ambiguous alternative gene
                     kch+=k*(nsame==4 ? 2 : 4);
                     newgene&=golg[nb[(nb1i>>(kch<<2))&0x7]];
                 }; return(newgene);
        case 6:  if ( selection<16) return(genegol[selection-8]);                        // choose gene with GoL encoding : needs fixing for selection 11,13
        case 7:  return(randnr);                                                         // choose random gene : should really update randnr outside to ensure indept
        default: fprintf(stderr,"Error in switch of ambiguous rotation resolution, should never reach here\n");
                 return(0ull);
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
    if(ancestor != rootgene)
        if((genedataptr = (genedata *) hashtable_find(&genetable, ancestor)) == NULL)
            fprintf(stderr,"error in hashaddgene, the ancestor %llx of gene %llx to be stored is not stored\n",ancestor,gene);
#endif
}
//.......................................................................................................................................................
extern inline void hashdeletegene(uint64_t gene,const char errorformat[]) {
#ifdef HASH
        if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) genedataptr->popcount--;
        else fprintf(stderr,errorformat,totsteps,gene);     // errorformat must contain %d and %llx format codes in this order
#endif
}
//.......................................................................................................................................................
extern inline void hashreplacegene(uint64_t gene1,uint64_t gene2,uint64_t ancestor,const char errorformat[]) {
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
    if(ancestor != rootgene)
        if((genedataptr = (genedata *) hashtable_find(&genetable, ancestor)) == NULL)
            fprintf(stderr,"error in hashreplacegene, the ancestor %llx of gene %llx to be stored is not stored\n",ancestor,gene2);
#endif
}
//.......................................................................................................................................................
extern inline void hashgeneextinction(uint64_t gene,const char errorformat[]) {
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
//.......................................................................................................................................................
extern inline void hashgeneactivity(uint64_t gene, const char errorformat[]) {
#ifdef HASH
        if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) genedataptr->activity ++;
        else fprintf(stderr,errorformat,4,totsteps,gene);
#endif
}
//------------------------------------------------------- hash quadtree inline fns ----------------------------------------------------------------------
extern inline uint64_t patt_hash(const uint64_t a, const uint64_t b, const uint64_t c, const uint64_t d) {
                                                        // this hash function works much better than that used in golly for example
    uint64_t a1,b1,c1,d1,r;
    a1 = a^0x7f0e1d2c3b4a5968ull;
    b1 = b^0xf0e1d2c3b4a59687ull;
    c1 = c^0xba9876543210fedcull;
    d1 = d^0x456789abcdef0123ull;
    r =  (a1>>13)+(b1>>23)+(c1>>29)+(d1>>31);
    r += ((d1<<16)+(c1<<8)+(b1<<4)+(a1<<2)) + (a1+b1+c1+d1);
    r += (r >> 11);
    // r=r*11400714819323198549ull;
    return r ;
}
//.......................................................................................................................................................
quadkey * newkey() {                                   // free keys kept in a linked list, allocated 1024 at a time
    quadkey *q;
    int i;
    uint64_t randnr;
    const int blocksize = 1024;
   
    if (freenodes == 0) {
        freenodes = (quadkey *) calloc(blocksize, sizeof(quadkey)) ;
        if (freenodes == 0) {
            fprintf(stderr,"No more memory can be allocated for quadkey chains.\n") ;
            exit(1);
        }
        quadchainmem += blocksize * sizeof(quadkey) ;
        for(i=0;i<(blocksize-1);i++) {
            freenodes[i].next = freenodes+i+1;
            RAND128P(randnr);
            freenodes[i].key = randnr;
        }
    }
    q = freenodes;
    freenodes = freenodes->next;
    return q ;
}
//.......................................................................................................................................................
extern inline uint64_t node_hash(const void *a, const void *b, const void *c, const void *d) {
   uint64_t r = (65537*(uint64_t)(d)+257*(uint64_t)(c)+17*(uint64_t)(b)+5*(uint64_t)(a));
   r += (r >> 11);
   return r ;
}
//.......................................................................................................................................................
extern inline quadnode * hash_patt_find(const uint64_t nw, const uint64_t ne, const uint64_t sw, const uint64_t se) {
#ifdef HASH
        quadnode *q;
        uint64_t h;
        int nr1,nr1s;
        h = patt_hash(nw,ne,sw,se);
        if((q = (quadnode *) hashtable_find(&quadtable, h)) != NULL) {
                                                        // check if pattern found in hash table is correct
            if( nw == q->nw.patt &&  ne == q->ne.patt && sw == q->sw.patt && se ==q->se.patt && !q->isnode) {
                q->hits++;q->active=1;
            }
            else {                                      // collision in hash table at leaf level
                quadcollisions++;
                fprintf(stderr,"at %d quadhash pattern collision %llx %llx %llx %llx hash %llx collides %llx %llx %llx %llx\n",
                        totsteps,nw,ne,sw,se,h,q->nw.patt,q->ne.patt,q->sw.patt,q->se.patt);
                while(q->next != NULL && (q = (quadnode *) hashtable_find(&quadtable, q->next->key)) != NULL) {
                // CODE DEVELOPMENT UP TO HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                }
            }
        }
        else {                                           // new node or pattern, save in hash table
            quadinit.isnode=0;
            quadinit.nw.patt=nw;
            quadinit.ne.patt=ne;
            quadinit.sw.patt=sw;
            quadinit.se.patt=se;
            quadinit.firsttime=totsteps;
            POPCOUNT64C(nw,nr1);nr1s=nr1;
            POPCOUNT64C(ne,nr1);nr1s+=nr1;
            POPCOUNT64C(sw,nr1);nr1s+=nr1;
            POPCOUNT64C(se,nr1);nr1s+=nr1;
            quadinit.pop1s =nr1s;
            hashtable_insert(&quadtable, h,(quadnode *) &quadinit);
            q = (quadnode *) hashtable_find(&quadtable, h);
        }
        return(q);
#else
        return(NULL);
#endif
}
//.......................................................................................................................................................
extern inline quadnode * hash_node_find(const quadnode * nw, const quadnode* ne, const quadnode* sw, const quadnode* se) {
#ifdef HASH
        quadnode *q;
        uint64_t h;
        h = node_hash(nw,ne,sw,se);
        if((q = (quadnode *) hashtable_find(&quadtable, h)) != NULL) {
            if(nw == q->nw.node && ne == q->ne.node && sw == q->sw.node && se == q->se.node && q->isnode) { // node found in hash table
                q->hits++;q->active=1;
                /* fprintf(stderr,"at %d quadhash node repeat %llx %llx %llx %llx hash %llx\n", totsteps,
                    (uint64_t) nw,(uint64_t) ne,(uint64_t) sw,(uint64_t) se,(uint64_t) h);  // simple recording for now, later do chaining or whatever */
            }
            else {                                       // collision in hash table at node level
                quadcollisions++;
                fprintf(stderr,"at %d quadhash node collision %llx %llx %llx %llx hash %llx collides %llx %llx %llx %llx\n", totsteps,
                    (uint64_t) nw,(uint64_t) ne,(uint64_t) sw,(uint64_t) se,(uint64_t) h,
                    (uint64_t) q->nw.node, (uint64_t) q->ne.node, (uint64_t) q->sw.node, (uint64_t) q->se.node);  // simple recording for now, later do chaining or whatever
            }
        }
        else {                                           // new node or pattern, save in hash table
            quadinit.isnode=1;
            quadinit.nw.node=(quadnode *)nw;
            quadinit.ne.node=(quadnode *)ne;
            quadinit.sw.node=(quadnode *)sw;
            quadinit.se.node=(quadnode *)se;
            quadinit.firsttime=totsteps;
            quadinit.pop1s=nw->pop1s+ne->pop1s+sw->pop1s+se->pop1s;
            hashtable_insert(&quadtable, h,(quadnode *) &quadinit);
            q = (quadnode *) hashtable_find(&quadtable, h);
            /* fprintf(stderr,"at %d quadhash node stored %llx %llx %llx %llx hash %llx\n", totsteps,
                    (uint64_t) nw,(uint64_t) ne,(uint64_t) sw,(uint64_t) se,(uint64_t) h);  // simple recording for now, later do chaining or whatever */
        }
        return(q);
#else
        return(NULL);
#endif
}
//.......................................................................................................................................................
quadnode * quadimage(uint64_t gol[],int N,int offset) {                                    // routine to generate a quadtree for an entire binary image of long words
                                                                                // makes use of global linear and quadratic size variable N2 for sizing internal arrays golp,golq
                                                                                // assumes that N is power of 2 and >=16
    unsigned int ij,ij1,n;
    uint64_t golp[N2>>6];
    quadnode * golq[N2>>8];
    extern void pack64neighbors(uint64_t gol[],uint64_t golp[],int N,int offset);
    
    pack64neighbors(gol,golp,N,offset);
    if(N<16) fprintf(stderr,"Error: attempting to quadtree image with size less than 16x16\n");
    n=N>>3;
                            /*  for(ij=0;ij<n*n;ij++) {
                                    if (ij%8 == 0) fprintf(stderr,"\n step %d ij %d",totsteps,ij);
                                    fprintf(stderr," %llx ",golp[ij]);
                                }
                                fprintf(stderr,"\n"); */
    for (ij=ij1=0;ij<n*n;ij+=2,ij1++) {
        golq[ij1]=hash_patt_find(golp[(ij+n)],golp[(ij+n)+1],golp[ij],golp[ij+1]); // hash_patt_find(nw,ne,sw,se) adds quad leaf entry if pattern not found
        ij+= ((ij+2)&(n-1)) ? 0 : n;                                            // skip odd rows since these are the northern parts of quads generated on even rows
    }
                            // fprintf(stderr,"8x8 patterns found\n");
    for(n >>= 1; n>1; n>>= 1) {
        for (ij=ij1=0;ij<n*n;ij+=2,ij1++) {
                            // fprintf(stderr,"n %d ij1 %d\n",n,ij1);
            golq[ij1]=hash_node_find(golq[(ij+n)],golq[(ij+n)+1],golq[ij],golq[ij+1]); // hash_node_find(nw,ne,sw,se) adds quad node entry if node not found
            ij+= ((ij+2)&(n-1)) ? 0 : n;                                        // skip odd rows since these are the northern parts of quads generated on even rows
        }
    }
    if(golq[0]!=NULL)
        if(golq[0]->hits > 1) fprintf(stderr,"image already found at t = %d\n",golq[0]->firsttime);
    return(golq[0]);
}
//.......................................................................................................................................................
void extract_components(short unsigned int label[],short unsigned int nlabel) {
    int i,ij;
    short unsigned int k,nside;

    for (ij=0;ij<N2;ij++) {
        if (label[ij] && !complist[label[ij]].first) {
            complist[label[ij]].first=ij+1;                                                             // ij+1 to avoid confusion of ij==0 with new component
            complist[label[ij]].E=N;
            complist[label[ij]].S=N;
        }
    }
    for (ij=0;ij<N2;ij++) {
        if (label[ij]) {
            if (complist[label[ij]].W > (ij&(N-1))) complist[label[ij]].W=(ij&(N-1));                   // adjust left boundary
            else if (complist[label[ij]].E < (ij&(N-1))) complist[label[ij]].E=(ij&(N-1));              // adjust right boundary
            if (complist[label[ij]].N > (ij>>log2N)) complist[label[ij]].N=(ij>>log2N);                 // adjust top boundary
            else if (complist[label[ij]].S < (ij>>log2N)) complist[label[ij]].S=(ij>>log2N);            // adjust bottom boundary
        }
    }
    for(i=1;i<nlabel;i++) {
        for(k=1;k<(1+complist[i].E-complist[i].W);k=k<<1);                                              // find nside as power of two square side length to fit component
        nside = k;                                                                                      // east-west minimum power of two side length
        for(k=1;k<(1+complist[i].S-complist[i].N);k=k<<1);
        nside = k > nside ? k : nside;                                                                  // larger of east-west and south-north minimum power of two side lengths
        complist[i].Nside = nside <16 ? 16 : nside;
        complist[i].corner = N*complist[i].N +complist[i].W;                                            // top left corner of component square image
        complist[i].quad = quadimage(gol,complist[i].Nside,complist[i].corner);                         // quadtree hash code of component image
    }
}
//------------------------------------------------------- pack012,3neighbors -------------------------------------------------------------------------------
#define deltaxy(ij,x,y)  (ij - (ij&Nmask) + (((ij+(x))&Nmask) + (y)*N)) & N2mask
//.......................................................................................................................................................
extern inline void pack012neighbors(uint64_t gol[],uint64_t golp[]) {              // routine to pack all up to 2nd neighbours in single word
    unsigned int ij,k;
    uint64_t gs;
    int nbx[24] = {-1,0,1,1,1,0,-1,-1,-2,-1,0,1,2,2,2,2,2,1,0,-1,-2,-2,-2,-2};
    int nby[24] = {-1,-1,-1,0,1,1,1,0,-2,-2,-2,-2,-2,-1,0,1,2,2,2,2,2,1,0,-1};

    for (ij=0;ij<N2;ij++)  golp[ij] = gol[ij];       // copy 1 bit gol to golp
    
    for (ij=0;ij<N2;ij++) {                          // copy up to 2nd neighbours to golp in upper 32-bit word
        for(k=0,gs=0ull;k<24;k++) {
            gs=gol[deltaxy(ij,nbx[k],nby[k])];
            golp[ij] |= gs<<(k+32);
        }
    }
}
//.......................................................................................................................................................
extern inline void pack0123neighbors(uint64_t gol[],uint64_t golp[]) {              // routine to pack all up to 3rd neighbours in single word
    unsigned int ij,k;
    uint64_t gs;
    int nbx[48] = {-1,0,1,1,1,0,-1,-1,-2,-1,0,1,2,2,2,2,2,1,0,-1,-2,-2,-2,-2,-3,-2,-1,0,1,2,3,3,3,3,3,3,3,2,1,0,-1,-2,-3,-3,-3,-3,-3,-3};
    int nby[48] = {-1,-1,-1,0,1,1,1,0,-2,-2,-2,-2,-2,-1,0,1,2,2,2,2,2,1,0,-1,-3,-3,-3,-3,-3,-3,-3,-2,-1,0,1,2,3,3,3,3,3,3,3,2,1,0,-1,-2};

    for (ij=0;ij<N2;ij++) golp[ij] = gol[ij];       // copy 1 bit gol to golp
    
    for (ij=0;ij<N2;ij++) {
        for(k=0,gs=0ull;k<48;k++) {
            gs=gol[deltaxy(ij,nbx[k],nby[k])];
            golp[ij] |= gs<<(k+8);                  // copy up to 2nd neighbours to golp in upper 7 bytes of word
        }
    }
}
//.......................................................................................................................................................
extern inline void pack49neighbors(uint64_t gol[],uint64_t golp[]) {              // routine to pack all up to 3rd neighbours in single word
    unsigned int ij,k;
    int nbx[6] = {1,0,2,0,-4,-4};
    int nby[6] = {0,1,0,2,0,-4};

    for (ij=0;ij<N2;ij++) golp[ij] = gol[ij];                                     // copy 1 bit gol to golp
    for(k=0;k<6;k++)                                                              // hierarchical bit copy and swap
        for (ij=0;ij<N2;ij++)
             golp[ij] |= gol[deltaxy(ij,nbx[k],nby[k])]<<(1<<k);                  // 8x8 packed arrays
    for (ij=0;ij<N2;ij++) golp[ij] = golp[ij]&0xfac8ffccfafaffffull;              // masks out 15 values in top row and left column to give 7x7 neighbourhoods
                                                                                  // mask removes bit numbers 16,18,24,26,32,33,36,37,48,49,50,52,53,56,58
}
//.......................................................................................................................................................
extern inline void pack16neighbors(uint64_t gol[],uint64_t golp[]) {              // routine to pack 4x4 subarrays into single words
    unsigned int ij,ij1,k;
    
    for (ij1=0;ij1<(N2>>6);ij1++) golp[ij1]=0ull;
    for(k=0;k<16;k++)
        for (ij=ij1=0;ij<N2;ij+=4,ij1++) {
             golp[ij1] |= gol[deltaxy(ij,k&0x3,k>>2)]<<k;                         // 4x4 packed arrays
             ij+= ((ij+4)&(N-1)) ? 0 : (4-1)*N;
        }
}
//.......................................................................................................................................................
extern void pack64neighbors(uint64_t gol[],uint64_t golp[],int n,int offset) {                     // routine to pack 8x8 subarrays into single words
    unsigned int ij,ij1,k,log2n;
    unsigned int n2 = n*n;
    for(log2n=k=1;k<n;k<<=1,log2n++);
    for (ij1=0;ij1<(n2>>6);ij1++) golp[ij1]=0ull;
    for(k=0;k<64;k++)
        for (ij=ij1=0;ij<n2;ij+=8,ij1++) {
             // golp[ij1] |= gol[deltaxy(ij,k&0x7,k>>3)]<<k;
             golp[ij1] |= gol[deltaxy(offset,(ij&(n-1))+(k&0x7),(ij>>log2n)+(k>>3))]<<k; // 8x8 packed arrays, assuming gol is lowest 1 bit only and golp length N2>>6
             ij+= ((ij+8)&(n-1)) ? 0 : (8-1)*n;
        }
}
//.......................................................................................................................................................
extern inline void compare_neighbors(uint64_t a[],uint64_t b[], int dx, int dy) {  // routine to compare packed pack neighbours with shift, result in a
    unsigned int ij,ijs,scrollN;
    uint64_t bij;
   
    if (last_scrolled) scrollN = N;
    else scrollN = 0;
    
    for (ij=0;ij<N2;ij++) {
        ijs=(ij-scrollN)&N2mask;
        bij=b[deltaxy(ijs,dx,dy)];
        a[ij] = (a[ij]|bij) ? a[ij]^bij : rootgene;
    }
}
//.......................................................................................................................................................
extern inline void compare_all_neighbors(uint64_t a[],uint64_t b[]) {  // routine to compare packed pack neighbours with shift, result in a
    unsigned int ij,ijs,scrollN;
    unsigned int d, k;
    uint64_t aij,bijk;
    int nbx[8] = {0,1,0,-1,1,1,-1,-1};   // N E S W NE SE SW NW
    int nby[8] = {-1,0,1,0,-1,1,1,-1};
    
    if (last_scrolled) scrollN = N;
    else scrollN = 0;
    
    for (ij=0;ij<N2;ij++) {
        ijs=(ij-scrollN)&N2mask;
        aij = a[ij];
        for (a[ij]=0ull,k=0;k<8;k++) {
            bijk=b[deltaxy(ijs,nbx[k],nby[k])];
            POPCOUNT64C((aij^bijk),d);
            d = (aij&&bijk) ? d : 0xff;
            a[ij]|=((uint64_t) d)<<(k<<3);
        }
    }
}
//.......................................................................................................................................................
extern inline void packandcompare(uint64_t newgol[],uint64_t working[],uint64_t golmix[]) {
    if (colorfunction==8) {
        pack49neighbors(newgol,working);
        if(offdx==0 && offdy==0 && offdt==0) {
            pack49neighbors(gol,golmix);
            compare_all_neighbors(golmix,working);  // compare all 8 directions N E S W NE SE SW NW
        }
        else {
            if (offdt<=-maxPlane) offdt=-maxPlane;
            if(offdt>0) offdt = 0;
            pack49neighbors(planesg[(newPlane-offdt)%maxPlane],golmix);
            compare_neighbors(golmix,working,offdx,offdy);                 // compare with a single direction (north) for gliders
        }
    }
}
//------------------------------------------------------- connected component labelling ---------------------------------------------------------------
// Combine two trees containing node i and j. Union by rank with halving - see https://en.wikipedia.org/wiki/Disjoint-set_data_structure
// Return the root of union tree in equivalence relation's disjoint union-set forest
short unsigned int lab_union(equivrec equiv[], short unsigned int i, short unsigned int j) {
    short unsigned int root = i;
    while (equiv[root].parent<root) {
        equiv[root].parent = equiv[equiv[root].parent].parent;        // halving algorithm to compress path on the fly
        root = equiv[root].parent;
    }
    if (i != j) {
        short unsigned int rooti = root;
        short unsigned int rootj;
        root = j;
        while (equiv[root].parent<root) {
            equiv[root].parent = equiv[equiv[root].parent].parent;   // halving algorithm to compress paths on the fly
            root = equiv[root].parent;
        }
        if (root == rooti) return root;
        rootj=root;
        if (equiv[rooti].rank < equiv[rootj].rank) {                 // swap rooti, rootj so that rooti has larger rank
            rootj=rooti;
            rooti=root;
        }
        equiv[rootj].parent = rooti;                                 // merge rootj into rooti
        if (equiv[rooti].rank == equiv[rootj].rank) {
            equiv[rooti].rank++;
        }
    }
    return root;
}
//.......................................................................................................................................................
extern inline short unsigned int label_cell(int ij,short unsigned int *nlabel) {
    short unsigned int clabel,clabel1,labelij;
    clabel=label[deltaxy(ij,0,-1)];                  // deltaxy takes periodic BCs into account (b in Suzuki Decision Tree, see Wu et.al.)
    if(clabel) labelij=equiv[clabel].parent;         // b
    else {                                           // N (=b) unlabelled
        clabel=label[deltaxy(ij,1,-1)];              // (c in Suzuki DT)
        if(clabel) {
            clabel1=label[deltaxy(ij,-1,-1)];        // NW (a in Suzuki DT)
            if(clabel1) {
                labelij=lab_union(equiv,clabel,clabel1); // c,a
            }
            else {
                clabel1=label[deltaxy(ij,-1,0)];     // W (d in Suzuki DT)
                if(clabel1) {
                    labelij=lab_union(equiv,clabel,clabel1); // c,d
                }
                else {
                    labelij=equiv[clabel].parent;    // c
                }
            }
        }
        else {
            clabel=label[deltaxy(ij,-1,-1)];         // NW (a in Suzuki DT)
            if(clabel) {
                labelij=equiv[clabel].parent;        // a
            }
            else {
                clabel=label[deltaxy(ij,-1,0)];      // W (d in Suzuki DT)
                if(clabel) {
                    labelij=equiv[clabel].parent;    // d
                }
                else {
                    clabel=++*nlabel;                // new label if a,b,c,d all without labels
                    equiv[clabel].parent=clabel;
                    labelij=clabel;
                    equiv[clabel].rank=0;
                }
            }
        }
    }
    return labelij;
}
//.......................................................................................................................................................
// Flatten the Union-Find tree and relabel the components.
void flattenlabels(equivrec equiv[],unsigned short int *nlabel) {
    short unsigned int xlabel = 0;
    short unsigned int i;
    for (i = 1; i < *nlabel + 1; i++)
    if (equiv[i].parent < i) {
        equiv[i].parent = equiv[equiv[i].parent].parent;
    }
    else {
        equiv[i].parent = xlabel;
        xlabel++;
    }
}
//.......................................................................................................................................................
// fast component labelling, updating global equivalence table equiv and placing labels in global array label
void label_components(uint64_t gol[]) {
    int ij;                                                 //   abc   scan mask to calculate label at position e
    short unsigned int nlabel = 0;                          //   de
    for(ij=0;ij<N2;ij++) label[ij]=0;
    for(ij=0;ij<(N2>>2);ij++) equiv[ij].parent=0;
    
    for(ij=0;ij<N2;ij++) {                                  // do first pass main array
        if (gol[ij]) {
            label[ij]=label_cell(ij,&nlabel);
        }
    }
    for(ij=0;ij<N;ij++) {                                   // redo first row sites to take periodic wrap around into account
        if (gol[ij]) {
            label[ij]=label_cell(ij,&nlabel);
        }
    }
    for(ij=N;ij<N2;ij+=N) {                                 // redo first column sites to take periodic wrap around into account
        if (gol[ij]) {
            label[ij]=label_cell(ij,&nlabel);
        }
    }
    flattenlabels(equiv,&nlabel);                           // single pass flatten of equivalence table
    for(ij=0;ij<N2;ij++) {                                  // do second pass of main array
        if (gol[ij]) {
            label[ij]=equiv[label[ij]].parent;
        }
    }
}
//------------------------------------------------------- geography -----------------------------------------------------------------------------------
void v_scroll(uint64_t newgol[],uint64_t newgolg[]) {
    int ij,scroll_needed;

    scroll_needed = 0;
    for (ij=N2-N;ij<N2;ij++) {                              // clear top row
        if(newgol[ij]) {
            scroll_needed = 1;
            break;
        }
    }
    last_scrolled = scroll_needed;                          // global variable needed for correct glider detection if scrolling
  
    for (ij=0;ij<N;ij++) {                                  // delete genes in bottom buffer row
        if(newgol[ij]) {
            newgol[ij]=0ull;
            hashdeletegene(newgolg[ij],"error in v_scroll hashdeletegene call for step %d with gene %llx\n");
            newgolg[ij]=gene0;
        }
    }
    
    if(!scroll_needed) return;
    
    for (ij=0;ij<N2-N;ij++) {                               // scroll rows down 1 leaving top row intact
        newgol[ij]=newgol[ij+N];
        newgolg[ij]=newgolg[ij+N];
    }
    for (ij=0;ij<N;ij++) {                                  // delete all states and genes in new bottom buffer row
        if(newgol[ij]) {
            newgol[ij]=0ull;
            hashdeletegene(newgolg[ij],"error in v_scroll hashdeletegene call for step %d with gene %llx\n");
            newgolg[ij]=gene0;
        }
    }
    for (ij=N2-N;ij<N2;ij++) {                              // clear top row
        if(newgol[ij]) {
            newgol[ij]=0ull;
            hashdeletegene(newgolg[ij],"error in v_scroll hashdeletegene call for step %d with gene %llx\n");
            newgolg[ij]=gene0;
        }
    }
}
//.......................................................................................................................................................
void random_soup(uint64_t gol[],uint64_t golg[],uint64_t newgol[],uint64_t newgolg[]) {
    int Nf,i,j,ij,i0,j0,i1,j1,d,k;
    uint64_t randnr,mask;
    static unsigned int rmask = (1 << 15) - 1;
    unsigned int density;
    
    if(rbackground) {                                           // homogeneous random background at rate rbackground/32768 per site per frame
        for(ij=0;ij<N2;ij++) {
            if((rand()&rmask) < rbackground) {
                if (!newgol[ij]) {                               // if live cell or multiplane, fill with game of life genome or random genome
                    if(randomsoup!=2 || selection < 8) {            // random gene if randomsoup non zero
                        RAND128P(randnr);
                        newgolg[ij] = randnr;
                        newgol[ij] = 1ull;
                    }
                    else {
                        newgolg[ij] = genegol[selection-8];         // otherwise GoL gene
                        newgol[ij] = 1ull;
                    }
                    hashaddgene(newgolg[ij],rootgene);
                }
            }
        }
        return;
    }
    
    if(randomsoup==2)
        if (totsteps & 0xf) return;                             // only execute once every 16 time steps
    
    density = initial1density;
    mask=(NbG==64 ? ~0 :(1ull<<NbG)-1ull);
    Nf = initfield;
    if (Nf==0 || Nf>N) Nf=N;
    i0 = j0 = (N>>1)-(Nf>>1);
    
    for (i=0; i<Nf; i++) {
        for (j=0; j<Nf; j++) {
            ij=i0+i+N*(j0+j);
            if(randomsoup==2) {                                 // border feathering as well as intermittent every 16 steps
                i1 = i<j ? i : j;                               // swap so that i1<=j1
                j1 = i<j ? j : i;
                d= j1< (Nf>>1) ? i1 : (i1 < Nf-j1 ? i1 : Nf-j1);// find Manhatten distance to border ij1
                density = (d <= 8 ? (initial1density >> (16-(d<<1))) : initial1density);
            }
            if(!newgol[ij]) {                                   // check whether this OK for selection 16-19 with genes everywhere and hashdelete
                if (selection<14) newgol[ij] = ((rand() & rmask) < density)?1ull:0ull;
                else if(selection>=16 && selection<=19) for (k=0;k<NbP;k++) newgol[ij] |= ((rand() & rmask) < density)?(1ull<<(k<<2)):0ull;
                else for (k=0;k<NbP;k++) newgol[ij] |= ((rand() & rmask) < density)?(1ull<<k):0ull;
                
                if (newgol[ij]) {                               // if live cell or multiplane, fill with game of life genome or random genome
                    if (selection<8) {
                        RAND128P(randnr);
                        newgolg[ij] = gene0^(randnr&mask);
                    }
                    else {
                        newgolg[ij] = genegol[selection-8];
                    }
                    hashaddgene(newgolg[ij],rootgene);
                }
            }
        }
    }
}
//------------------------------------------------------- update -----------------------------------------------------------------------------------
void update(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]){
    /* update GoL for toroidal field which has side length which is a binary power of 2 */
    /* encode without if structures for optimal vector treatment */
    int s, s0, k, k1, kodd, nmut, crot;
    unsigned int kch,rulemodij,add2nd,mask1st;
    unsigned int select23live,pos_canon_neutral,survival,overwrite,enforcebirth,add2ndmask1st,nongolnottwice;
    int nb[8], nbc, nbch, ij, i, j, jp1, jm1, ip1, im1;
    uint64_t g, gs, nb1i, nb2i, randnr, r2;
    uint64_t nbmask, nbmaskr;
    uint64_t newgene, ancestor, livegenes[3];
    uint64_t s2or3, birth, statflag, nextgolstate;

    survival = survivalmask;
    overwrite = overwritemask;
    select23live = ((repscheme & R_0_2sel_3live)?1:0)+((repscheme & R_1_2sel_2live)?2:0);
    pos_canon_neutral = ((repscheme & R_2_canonical_nb)?1:0)+((repscheme & R_3_neutral_pos)?2:0);
    enforcebirth =((repscheme & R_4_enforce3birth)?1:0)+((repscheme & R_5_enforce2birth)?2:0);
    add2ndmask1st = ((repscheme & R_6_2ndnb_genes)?1:0)+((repscheme & R_7_1stnb_masks)?2:0);
    nongolnottwice = ((repscheme & R_8_nongolstat)?1:0)+((repscheme & R_9_nongolstatnbs)?2:0);
    add2nd = add2ndmask1st&0x1;

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
        rulemodij = (rulemod&0x2) ? (ij>=(N2>>1) ? 1 : 0) : (rulemod&0x1);   // if rulemod bit 1 is on then split into half planes with/without mod
        nbmask = 0;
        if(s>1) {
          if (repscheme & R_17_quadrant_2nb1)  add2ndmask1st =     (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
          if (repscheme & R_18_quadrant_ngol)  nongolnottwice =    (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
          add2nd = add2ndmask1st&0x1; mask1st=(add2ndmask1st>>1)&0x1;
          
          if(nongolnottwice&0x1) {                                          // check for non GoL changed states
            if((nongolnottwice>>1)&0x1) {                                   // check in neighborhood
                int sng;                                                    // sum of nbs with non GoL rule change bit set
                for (k=0,sng=0;k<8;k++) {                                   // calc number of neighbours resulting from nongol rule change
                    sng += (golgstats[nb[k]]&F_nongolchg)?1:0;              // neighbors with state set by a non GoL change
                }
                if(sng) rulemodij = 0ull;
            }
            if(golgstats[ij]&F_nongolchg) rulemodij = 0ull;                 // if central state the result of a non GoL rule change
          }
          if(mask1st&&rulemodij) {                                         // recalculate effective new s as less than original sum sm
            for (k=0,nbmaskr=0;k<8;k++) {                                   // depending on rotated overlay of masks in live neighbour genes
                if(gol[nb[k]]) {
                    g= golg[nb[k]];                                         // fnal gene has ncoding 0s then 8 bit mask
                    //if(!((g>>8) & codingmask)) nbmaskr |= g&0xff;           // tried |=, &=, ^= .
                    if(!((g>>8) & codingmask)) nbmaskr |= (0x2<<(g&0x7L))&0xff;  // only one bit on for gene masks in this version: not self, so 0x2L

                }
                nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);           // 8 bit rotate right
            }
            for (k=0,s=0,nb1i=0ull;k<8;k++) {                               // recalculate sum and nb1i using combined mask
                gs =gol[nb[k]] & (((~nbmaskr)>>k)&0x1ull);                  // if mask bit set, count as if dead
                nbmask |= gs<<k;                                            // also calculate nbmask for use below
                s += gs;
                nb1i = (nb1i << (gs<<2)) + (gs*k);
            }
            if(s!=s0) s2or3 = (s>>2) ? 0ull : (s>>1);                       // redo s2or3 calculation for modified s
          }
          else for (k=0,nbmask=0;k<8;k++) nbmask |= (gol[nb[k]]<<k);        // 8-bit mask of GoL states of 8 nbs, clockwise from top left
        } // end if s>1
        if (s2or3) {                                                        // if 2 or 3 neighbours alive
            if (repscheme & R_quadrant) {                                   // quarter the plane with 4 different parameter values
                if (repscheme & R_14_quadrant_sele)  select23live =      (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
                if (repscheme & R_15_quadrant_posn)  pos_canon_neutral = (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
                if (repscheme & R_16_quadrant_enfb)  enforcebirth =      (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
                if (repscheme & R_19_quadrant_surv)  survival =          (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
                if (repscheme & R_20_quadrant_over)  overwrite =         (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
                canonical = pos_canon_neutral&0x1;       // global value since needed in ...difft2-6 subroutines
            }
            birth = 0ull;
            newgene = 0ull;

            if (s&0x1ull) {  // s==3                                        // allow birth (with possible overwrite)
              statflag |= F_3_live;                                         // record instance of 3 live nbs
              if ((0x1ull&overwrite)|(0x1ull&~gol[ij]) ) {                  // central site empty or overwrite mode
                birth = 1ull;                                               // birth flag
                for(k=0;k<s;k++) livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]]; // live neighbour genes
                statflag |= F_3_livenbs & (nbmask<<16);                     // record live neighbour pattern
                if((livegenes[0]^livegenes[1])|(livegenes[0]^livegenes[2])) { // genes not all same, need ancestor calculation
                  kch=selectdifft3(nbmask, &crot, &kodd);nbch=nb[kch];
                  if (select23live&0x1) {                                   // execute selective replication of one of two otherwise unchosen live genes
                      nb2i = 0ull;
                      for(k1=k=0;k<3;k++) {                                 // choice of two other live genes for possible ancestor
                          nbc=(nb1i>>(k<<2))&0x7;
                          if(nb[nbc]!=nbch) nb2i = (nbc<<(k1++<<2))+nb2i;
                      }
                      if (add2nd) selectone_nbs(s,nb2i,nb,gol,golg,&birth,&newgene);  //2nd nb modulation
                      else selectone_of_2(s,nb2i,nb,golg,&birth,&newgene,kch);
                      if (birth==0ull) {                                    // optional reset of ancestor & birth if no ancestors chosen in selectone
                        if((enforcebirth&0x1)||rulemodij)  {                // birth cannot fail or genes don't matter or no modification to gol rules
                            newgene = golg[nbch];
                            birth = 1ull;
                        }
                      }
                      else statflag |= F_2select;                           // ancestor has been chosen in selectone_of_2
                  }
                  else {
                      newgene = golg[nbch];
                  }
                } // end if not all live neighbors the same
                else {
                    statflag |= F_3g_same;
                    newgene = livegenes[0];                                 // genes all the same : copy first one
                    if((~enforcebirth&0x1) && rulemodij) birth = 0ull;      // no birth for 3 identical genes if not enforcebirth3 and rulemod
                }
              } // end central site empty or overwrite mode
            }  // end if s==3
            else {  // s==2                                                 // possible birth as exception to GoL rule
                statflag |= F_2_live;
                if (((select23live>>1)&0x1)&&(rulemodij||gol[ij])) {       // rule departure from GOL allowed or possible overwrite
                    if ((0x1ull&(overwrite>>1))||!gol[ij]) {                // either overwrite on for s==2 or central site is empty
                        if (add2nd) selectone_nbs(s,nb1i,nb,gol,golg,&birth,&newgene); //2nd nb modulation
                        else {
                            nbmask = (0x1ull<<(nb1i&0x7)) + (0x1ull<<((nb1i>>4)&0x7));
                            kch=selectdifft2(nbmask, &crot, &kodd);
                            if (repscheme & ((pos_canon_neutral>>1)&0x1)) { // enforce gene independent birth
                                newgene = golg[nb[kch]];
                                birth = 1ull;
                            }
                            else selectone_of_2(s,nb1i,nb,golg,&birth,&newgene,kch);
                            if(repscheme & R_10_13_2birth_k4) {             // birth only on the active subset of 4 canonical 2-live nb configs if any are active
                                if(~repscheme & (R_10_2birth_k0<<(kch-1))) birth = 0ull;   // cancel birth if chosen configuration not active
                            }
                        }

                        if(!birth && (enforcebirth&0x2)) {
                            nbmask = (0x1ull<<(nb1i&0x7)) + (0x1ull<<((nb1i>>4)&0x7));
                            kch=selectdifft2(nbmask, &crot, &kodd);
                            newgene = golg[nb[kch]];
                            if(repscheme & R_10_13_2birth_k4) {
                                if(repscheme & (R_10_2birth_k0<<(kch-1))) birth = 1ull;
                                else birth = 0ull;
                            }
                            else birth = 1ull;
                        }
                        if (birth) statflag |= F_2select;
                    }
                }
            }

            if(birth){
                // compute random events for multiple single bit mutations, as well as mutation position nmut
                r2=1ull;
                ancestor = newgene;
                while (r2) {
                    RAND128P(randnr);                                           // inline exp so compiler recognizes auto-vec,
                    r2 = randprob(pmutmask,(unsigned int) randnr);
                    nmut = (randnr >> 56) & 0x3f;                               // choose mutation position for length 64 gene : from bits 56:61 of randnr
                    newgene = newgene ^ (r2<<nmut);                             // introduce single mutation with probability pmut = probmut
                    statflag = statflag | F_mutation;
                }

                if(gol[ij]) {                                               // central old gene present: overwritten
                    hashdeletegene(golg[ij],"step %d hash delete error 1 in update, gene %llx not stored\n");
                }
                hashaddgene(newgene,ancestor);

                newgol[ij]  =  1ull;                                        // new game of life cell value: alive
                newgolg[ij] =  newgene;                                     // if birth then newgene
                statflag = statflag | F_birth;
            } // end birth
            else {
                if ((survival&s&0x1ull)|((survival>>1)&(~s)&0x1ull)|((~rulemodij)&0x1ull)) {// (surv bit 0 and s==3) or (surv bit 1 and s==2) or not rulemod1ij
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

    if(randomsoup) random_soup(gol,golg,newgol,newgolg);
    if(vscrolling) v_scroll(newgol,newgolg);
    if (colorfunction==8) packandcompare(newgol,working,golmix);
    else if (colorfunction==9) label_components(newgol);
    
    for (ij=0; ij<N2; ij++) {       // complete missing hash table records of extinction and activities
        if(gol[ij]) hashgeneextinction(golg[ij],"hash extinction storage error in update, gene %llx not stored\n");
        if(newgol[ij]) hashgeneactivity(newgolg[ij],"hash activity storage error in update, gene %llx not stored\n");
    }
    if(quadtree) qimage = quadimage(newgol,N,0);

}
//------------------------------------------------------- update_lut_sum -----------------------------------------------------------------------------------
void update_lut_sum(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]){    // selection models 8,9
// update GoL for toroidal field which has side length which is a binary power of 2
// encode without if structures for optimal vector treatment
/*
        0 <= s <= 8
        genome =
        8 bits for each possible s value for birth (center site 0) : exclude s=0 no birth of isolated
        8 bits for each possible s value for survival/death (center site 1) : exclude s=0 no survival of isolated
        (assume s = 0 => 0 always)
        using b0 for s=0, b1 for s=1, etc
              76543210    76543210
        GOL = 00000100(0)|00000110(0)
            = 0000 0100 0000 0110
            = 0x0406
                                                                                                                */
    int s, s2or3, k, nmut, kch, kodd, crot, nbest, nsame;
    unsigned int rulemodij;
    int nb[8],  ij, i, j, jp1, jm1, ip1, im1;
    uint64_t genecode, gols, nb1i, nbmask, found, randnr, r2;
    uint64_t newgene, ancestor;
    uint64_t  survive, birth, overwrite,survivalgene,smask, bmask, statflag, ncodingmask, allcoding;

    canonical = repscheme & R_2_canonical_nb;                                   // set global choice of canonical rotation bit choice for selectdifftx
    survivalgene = repscheme & R_0_survivalgene;                                // gene determining survival is 1: central gene 2: determined by neighbours
    smask = (uint64_t) survivalmask;                                            // convert to 64 bit mask for efficient usage here
    bmask = (uint64_t) birthmask;
    ncodingmask = (1ull<<ncoding)-1ull;                                         // mask for number of bits coding for each lut rule: <=4 for birth and survival case
    if (ncoding==4) allcoding = 0xffffffffffffffff;                             // mask for total gene coding region
    else if (ncoding ==2) allcoding = 0xffffffff;
    else allcoding = 0xffff;
    for (ij=0; ij<N2; ij++) {                                                   // loop over all sites of 2D torus with side length N
        i = ij & Nmask;  j = ij >> log2N;                                       // row & column
        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                       // toroidal (j+1)*N and (j-1)*N
        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                             // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;               // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (s=0,nb1i=0ull,k=0;k<8;k++) {                                       // packs non-zero nb indices in first up to 8*4 bits
            gols=gol[nb[k]];                                                    // whether neighbor is alive
            s += gols;                                                          // s is number of live nbs
            nb1i = (nb1i << (gols<<2)) + (gols*k);                              // nb1i is packed list of live neighbour indices
        }
        statflag = 0ull;
        birth = 0ull;
        if (s) {
            s2or3 = (s>>2) ? 0ull : (s>>1);                                     // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
            gols = s2or3 ? (gol[ij] ? 1ull : (s&1ull ? 1ull : 0ull )) : 0ull;   // GoL calculation next state for non-genetic gol plane
            rulemodij = (rulemod&0x2) ? (ij>=(N2>>1) ? 1 : 0) : (rulemod&0x1);  // if rulemod bit 1 is on then split into half planes with/without mod
            if(rulemodij) {
                overwrite = overwritemask&(0x1ull<<(s-1));
                if (selection==9) {                                             // selection == 9
                    for (genecode=0ull,k=0;k<8;k++) {                           // decodes genes with variable length encoding
                        if (gol[nb[k]]) {
                            if (!survivalgene && gol[ij]) {
                                PATTERN4(golg[nb[k]], (s&0x7), found);          //survival?
                                genecode |= found? 1ull << (s-1) : 0ull;
                            }
                            if (overwrite || !gol[ij]) {
                                PATTERN4(golg[nb[k]], (s|0x8), found);          //birth?
                                genecode |= found? 1ull << (s-1+8) : 0ull;
                            }
                        }
                    }
                    if (survivalgene && gol[ij]) {                              // survival determined by central gene in this case
                                PATTERN4(golg[ij], (s&0x7), found);             // survival?
                                genecode |= found? 1ull << (s-1) : 0ull;
                    }
                    survive = ((genecode&smask)>>(s-1)) & 0x1ull;
                    genecode>>=8;
                    birth   = ((genecode&bmask)>>(s-1)) & 0x1ull;
                }
                else {                                                          // selection == 8
                    if(repscheme&R_1_nb_OR_AND)
                        for (genecode=0ull,k=0;k<8;k++)                         // decodes genes with fixed length encoding by OR
                            genecode |= (gol[nb[k]]?golg[nb[k]]:0ull);          // OR of live neighbours encodes birth rule & survival rule
                    else
                        for (genecode=allcoding,k=0;k<8;k++)                    // decodes genes with fixed length encoding by AND
                            genecode &= (gol[nb[k]]?golg[nb[k]]:allcoding);     // AND of live neighbours encodes birth rule & survival rule
                    if(survivalgene) genecode = (genecode&(0xffffffffull<<32)) | (golg[ij]&0xffffffffull);  // if central gene determines survival
                    survive=(((genecode>>((s-1)*ncoding)) & ncodingmask) == ncodingmask) && ((smask>>(s-1))&1ull) ? 1 : 0;
                    if (overwrite || !gol[ij]) birth=(((genecode>>((8+(s-1))*ncoding)) & ncodingmask) == ncodingmask) && ((bmask>>(s-1))&1ull) ? 1 : 0;
                    // if(s ==3 && genecode!=genegol[selection-8]) fprintf(stderr,"genecode %llx != genegol %llx at ij %d\n",genecode,genegol[selection-8],ij);
                }
            }
            else {
                    survive = s2or3;
                    birth = s2or3&s&0x1ull&~gol[ij];
            }
        }
        else survive = birth = s2or3 = gols = 0ull;

        if(birth) { // birth
            RAND128P(randnr);                                               // inline exp so compiler recognizes auto-vec,
            if (repscheme & R_7_random_resln) {
                kch = ((randnr>>32)&0xffff) % s;                            // choose random in this option only
                newgene = golg[nb[(nb1i>>(kch<<2))&0x7]];
            }
            else {
                nbest=selectone_of_s(s,nb1i,nb,golg,&birth,&newgene,&nbmask);// selection scheme depends on repscheme parameter
                if(nbest>1) {
                    kch=selectdifft(nbest,nbmask,&crot,&kodd,&nsame);        // kch is chosen nb in range 0-7, nsame gives the number of undistinguished positions in canonical rotation
                    if(nsame) newgene = disambiguate(kch, nb1i, nb, golg, nsame, &birth, randnr); // restore symmetry via one of 8 repscheme options
                    else newgene = golg[nb[kch]];
                    // fprintf(stderr,"lut step %d ij %d s %d nbest %d nbmask %llx nsame %d crot %d newgene %llx nb1i %llx kch %d gol[nb[kch]] %llx\n",totsteps, ij, s, nbest, nbmask, nsame, crot, newgene, nb1i, kch, gol[nb[kch]]);
                }
            }
            if (birth) {                                                    // ask again because disambiguate may turn off birth
                statflag |= F_birth;
                r2=1ull;                                                    // compute random events for single bit mutation, as well as mutation position nmut
                ancestor = newgene;
                while (r2) {
                    RAND128P(randnr);                                       // inline exp so compiler recognizes auto-vec,
                    r2 = randprob(pmutmask,(unsigned int) randnr);
                    nmut = (randnr >> 56) & 0x3f;                           // choose mutation position for length 64 gene : from bits 56:61 of randnr
                    newgene = newgene ^ (r2<<nmut);                         // introduce single mutation with probability pmut = probmut
                    statflag = statflag | F_mutation;
                }
                newgol[ij]  =  1ull;                                        // new game of life cell value: alive
                newgolg[ij] =  newgene;                                     // if birth then newgene
                if(gol[ij]) hashdeletegene(golg[ij],"step %d hash delete error 1 in update_lut_sum, gene %llx not stored\n");
                hashaddgene(newgene,ancestor);
            }
            else {                                                          // no birth, stays as before
                newgol[ij] = gol[ij];
                newgolg[ij] = golg[ij];
                if(gol[ij] && survive) statflag |= F_survival;              // could be failed overwrite birth attempt and survival
            }
        }
        else if(gol[ij]) {                                                       // death/survival
            if(survive) {                                                   // survival coded
                statflag |= F_survival;
                newgol[ij]  = gol[ij];                                      // new game of life cell value same as old
                newgolg[ij] = golg[ij];                                     // gene stays same
            }
            else {                                                          // death
                statflag |= F_death;
                newgol[ij]  = 0ull;                                         // new game of life cell value dead
                newgolg[ij] = 0ull;                                         // gene dies
                hashdeletegene(golg[ij],"step %d hash delete error 2 in update, gene %llx not stored\n");
            }
        }
        else {                                                              // empty and no birth, stays empty
            newgol[ij] = gol[ij];
            newgolg[ij] = golg[ij];
        }

        if(newgol[ij]!=gols) {
            statflag |= F_notgolrul;
            if(newgol[ij]) statflag |= F_nongolchg;
        }
        if(gol[ij]) statflag |= F_golstate;                                 // this is the last gol state, not the new state
        newgolgstats[ij] = statflag;
    }  // end for ij

    if(randomsoup) random_soup(gol,golg,newgol,newgolg);
    if(vscrolling) v_scroll(newgol,newgolg);
    if (colorfunction==8) packandcompare(newgol,working,golmix);

    for (ij=0; ij<N2; ij++) {       // complete missing hash table records of extinction and activities
        if(gol[ij]) hashgeneextinction(golg[ij],"hash extinction storage error in update, gene %llx not stored\n");
        if(newgol[ij]) hashgeneactivity(newgolg[ij],"hash activity storage error in update, gene %llx not stored\n");
    }
    if(quadtree) qimage = quadimage(newgol,N,0);
}
//------------------------------------------------------------- update_lut_dist -----------------------------------------------------------------------------------
void update_lut_dist(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]){     // selection models 10,11
// update GoL for toroidal field which has side length which is a binary power of 2
// encode without if structures for optimal vector treatment
/*
    configurations are distinguished by number of ones in two classes (corner,edge-centred) of the 8-bit live neighbour pattern
    for the different s values                 0  1  2  3  4  5  6  7  8
    there are a nr of partitions               1  2  3  4  5  4  3  2  1   total 25   with se the number in edge centred sites NSEW
    chosen nr of partitions for exploration    0  2  3  4  5  4  3  2  0   total 23
    if we exclude the cases s = 0,8, then there are 23 bits required to distinguish these cases
                                                                                                                */
    int s, s1, se, s0, s2or3, k, kodd, nmut, crot, found, nbest, nsame;
    uint64_t survive,birth,overwrite,survivalgene,smask,bmask,rulemodij;
    static uint64_t summasks[7] = {0x3ull,0x7ull,0xfull,0x1full,0xfull,0x7ull,0x3ull};
    static int sumoffs[7] = {0,2,5,9,14,18,21};
    int nb[8], ij, i, j, jp1, jm1, ip1, im1;
    unsigned int kch=0;
    uint64_t genecode, gols, nb1i, nbmask, randnr, r2;
    uint64_t gene, newgene, ancestor;
    uint64_t statflag;
    
    canonical = repscheme & R_2_canonical_nb;                                      // set global choice of canonical rotation bit choice for selectdifftx
    survivalgene = repscheme & R_0_survivalgene;                                   // gene determining survival is 1: central gene 2: determined by neighbours
    smask = (uint64_t) survivalmask;                                               // convert to 64 bit mask for efficient usage here
    bmask = (uint64_t) birthmask;

    for (ij=0; ij<N2; ij++) {                                                      // loop over all sites of 2D torus with side length N
        i = ij & Nmask;  j = ij >> log2N;                                          // row & column
        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                          // toroidal (j+1)*N and (j-1)*N
        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                                // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;                  // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (s=se=0,nb1i=0ull,k=0;k<8;k++) {                                       // packs non-zero nb indices in first up to 8*4 bits
            gols=gol[nb[k]];                                                       // whether neighbor is alive
            s += gols;                                                             // s is number of live nbs
            se += k&0x1&gols;                                                      // se is number of edge-centred live neighbours (odd k)
            nb1i = (nb1i << (gols<<2)) + (gols*k);                                 // nb1i is packed list of live neighbour indices
        }

        statflag = newgene = survive = birth = 0ull;
        if (s>0 && s<8) {
            s1 = s-1;
            s2or3 = (s>>2) ? 0ull : (s>>1);                                        // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
            gols = s2or3 ? (gol[ij] ? 1ull : (s&1ull ? 1ull : 0ull )) : 0ull;      // GoL calculation next state for non-genetic gol plane
            rulemodij = (rulemod&0x2) ? (ij>=(N2>>1) ? 1 : 0) : (rulemod&0x1);     // if rulemod bit 1 is on then split into half planes with/without mod
            if (rulemodij) {                // NB need to put gene calculation outside that we cna do genetic propagation with GoL rulemod off
                overwrite = overwritemask&(0x1ull<<s1);
                overwrite = (overwrite || !gol[ij]) ? 1ull : 0ull;
                if (gol[ij]) survive = (smask>>sumoffs[s1])&summasks[s1] ? 1ull : 0ull;
                if (overwrite) birth = (bmask>>sumoffs[s1])&summasks[s1] ? 1ull : 0ull;
                if (survive|birth) {
                    s0 = (s-4 > 0 ? s-4 : 0);
                    if (gol[ij])   survive = (smask>>(sumoffs[s1]+se-s0))&0x1ull;   // refine decisions for specific combination of s and se
                    if (overwrite) birth = (bmask>>(sumoffs[s1]+se-s0))&0x1ull;     // only allowed if birth,survivalmask permits (ask this before consulting genes)
                    if (survive|birth) {                                            // complete determination of birth or survival
                        if (selection==11) {
                            survive=birth=0ull;
                            for (k=0;k<8;k++)  {                                    // decodes genes with variable length encoding only for current s,se
                                if (gol[nb[k]]) {                                   // combine information from genes of all live neighbours
                                    gene = golg[nb[k]] & (0xf0f0f0f0f0f0f0f0ull | (0x0101010101010101ull<<(se-s0))); // focus gene down to specific required se bit for exact matching
                                    if (!survivalgene && gol[ij]) {                 // coding is 8 bits [(b/s) (s1 2 1 0) (se subset mask 3 2 1 0)]  (exception case s==4,se==4 is x0000001)
                                        PATTERN8(gene, (s==4 && se==4) ? 0x01 : ((s<<4)|(1ull<<(se-s0))), found); //survival rule found? final decision for survival
                                        survive |= found? 1ull : 0ull;                                                   // special case codes for 5th case se==4 for s==4
                                    }
                                    if (overwrite) {
                                        PATTERN8(gene, (s==4ull && se==4ull) ? 0x81 :(((8|s)<<4)|(1ull<<(se-s0))), found); //birth rule found? final decision for birth
                                        birth |= found? 1ull : 0ull;
                                    }
                                }
                            }
                            if (survivalgene && gol[ij]) {                                                           // survival determined by central gene in this case
                                    PATTERN8(golg[ij], (s==4 && se==4) ? 0x01 : ((s<<4)|(1ull<<(se-s0))), found);    //survival rule found? final decision for survival
                                    survive |= found? 1ull : 0ull;                                                   // special case codes for 5th case se==4 for s==4
                            }
                        }
                        else {                                                      // selection == 10
                            if(repscheme&R_1_nb_OR_AND)
                                for (genecode=0ull,k=0;k<8;k++)                     // decodes genes with fixed length encoding by OR
                                    genecode |= (gol[nb[k]]?golg[nb[k]]:0ull);      // OR of live neighbours encodes birth rule & survival rule
                            else
                                for (genecode=~0ull,k=0;k<8;k++)                    // decodes genes with fixed length encoding by AND
                                    genecode &= (gol[nb[k]]?golg[nb[k]]:~0ull);     // AND of live neighbours encodes birth rule & survival rule
                            if(survivalgene) genecode = (genecode&(0xffffffffull<<32)) | (golg[ij]&0xffffffffull);  // if central gene determines survival
                            genecode &= (bmask << 32)|smask;                        // not required as now already tested above
                            if (gol[ij]) survive &= (genecode>>(sumoffs[s1]+(se-s0)))&0x1ull;
                            if (overwrite) birth &= (genecode>>(32+sumoffs[s1]+(se-s0)))&0x1ull;
                        }
                    }
                }
            }
            else {
                survive = s2or3;
                birth = s2or3&s&0x1&~gol[ij];
            }
        }
        else survive = birth = s2or3 = gols = 0ull;

        if(birth) {
            RAND128P(randnr);                                                   // inline exp so compiler recognizes auto-vec,
            if (repscheme & R_7_random_resln) {
                kch = ((randnr>>32)&0xffff) % s;                                // choose random in this option only
                newgene = golg[nb[(nb1i>>(kch<<2))&0x7]];
            }
            else {
                nbest=selectone_of_s(s, nb1i, nb, golg, &birth, &newgene, &nbmask); // selection scheme depends on repscheme parameter
                if(nbest>1) {
                    kch=selectdifft(nbest, nbmask, &crot, &kodd, &nsame);
                    if(nsame) newgene = disambiguate(kch, nb1i, nb, golg, nsame, &birth, randnr); // restore symmetry via one of 8 repscheme options
                    else newgene = golg[nb[kch]];
                }
            }
            if (birth) {                                                        // renewed query as possibly updated in disambiguate
                statflag |= F_birth;
                r2=1ull;                                                        // compute random events for single bit mutation, as well as mutation position nmut
                ancestor = newgene;
                while (r2) {                                                    // allow multiple point mutations
                    RAND128P(randnr);                                           // inline exp so compiler recognizes auto-vec,
                    r2 = randprob(pmutmask,(unsigned int) randnr);
                    nmut = (randnr >> 56) & 0x3f;                               // choose mutation position for length 64 gene : from bits 56:61 of randnr
                    newgene = newgene ^ (r2<<nmut);                             // introduce single mutation with probability pmut = probmut
                    statflag = statflag | F_mutation;
                }
                newgol[ij]  =  1ull;                                            // new game of life cell value: alive
                newgolg[ij] =  newgene;                                         // if birth then newgene
                if(gol[ij]) hashdeletegene(golg[ij],"step %d hash delete error 1 in update_lut_sum, gene %llx not stored\n");
                hashaddgene(newgene,ancestor);
            }
            else {                                                              // no birth, stays as before
                newgol[ij] = gol[ij];
                newgolg[ij] = golg[ij];
                if (gol[ij]) statflag |= F_survival;                            // in case of failed overwrite, survival is an option
            }
        }
        else if (gol[ij]) {                                                               // death/survival   (careful! genecode changed during this processing)
            if(survive) {
                statflag |= F_survival;
                newgol[ij]  = gol[ij];                                              // new game of life cell value same as old
                newgolg[ij] = golg[ij];                                             // gene stays same
            }
            else {                                                                  // gene bit = 0 => death
                statflag |= F_death;
                newgol[ij]  = 0ull;                                                 // new game of life cell value dead
                newgolg[ij] = 0ull;                                                 // gene dies
                hashdeletegene(golg[ij],"step %d hash delete error 2 in update, gene %llx not stored\n");
            }
        }
        else{                                                                       // stays empty
            newgol[ij] = gol[ij];
            newgolg[ij] = golg[ij];
        }
        
        if(newgol[ij]!=gols) {
            statflag |= F_notgolrul;
            if(newgol[ij]) statflag |= F_nongolchg;
        }

        if(gol[ij]) statflag |= F_golstate;                                         // this is the last gol state, not the new state
        newgolgstats[ij] = statflag;
    }  // end for ij

    if(randomsoup) random_soup(gol,golg,newgol,newgolg);
    if(vscrolling) v_scroll(newgol,newgolg);
    if (colorfunction==8) packandcompare(newgol,working,golmix);

    for (ij=0; ij<N2; ij++) {       // complete missing hash table records of extinction and activities
        if(gol[ij]) hashgeneextinction(golg[ij],"hash extinction storage error in update, gene %llx not stored\n");
        if(newgol[ij]) hashgeneactivity(newgolg[ij],"hash activity storage error in update, gene %llx not stored\n");
    }
    if(quadtree) qimage = quadimage(newgol,N,0);
}
//------------------------------------------------------- update_lut_canon_rot -----------------------------------------------------------------------------------
void update_lut_canon_rot(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]){     // selection models 12,13
/*
    configurations are distinguished by number of ones and the canonical rotation of the 8-bit live neighbour pattern
    the canonical rotation is the rotation with the minimum numerical value as 8-bit number (bits are numbered clockwise from 0-7)
    for the different s values                 0  1  2  3  4  5  6  7  8
    there are a nr of canonical rotations      1  1  4  7 10  7  4  1  1   total 36
    if we exclude the cases s=0,1,7,8 which can create growth artefacts, then there are 32 bits
                                                                                                                */
    int s, smid, s2, s2or3, k, k1, kodd, nmut, crot, crot5, crotmod5, pat, found, nbest, nsame;
    uint64_t survive,birth,overwrite,survivalgene,smask,bmask,rulemodij;
    static uint64_t summasks[5] = {0xfull,0x7full,0x3ffull,0x7full,0xfull};
    static int sumoffs[5] = {0,4,11,21,28};
    int nb[8], ij, i, j, jp1, jm1, ip1, im1;
    unsigned int kch=0;
    uint64_t genecode, genecode1, gols, nb1i, nbmask, randnr, r2;
    uint64_t newgene, ancestor;
    uint64_t statflag;
    
    canonical = repscheme & R_2_canonical_nb;                                  // set global choice of canonical rotation bit choice for selectdifftx
    survivalgene = repscheme & R_0_survivalgene ? 1ull : 0ull;                 // gene determining survival is 1: central gene 2: determined by neighbours
    smask = (uint64_t) survivalmask;                                           // 32 bits of survivalmask used to limit space of rules, convert to 64 bit masks for efficient usage here
    bmask = (uint64_t) birthmask;                                              // 32 bits of birthmask

    for (ij=0; ij<N2; ij++) {                                                  // loop over all sites of 2D torus with side length N
        i = ij & Nmask;  j = ij >> log2N;                                      // row & column
        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                      // toroidal (j+1)*N and (j-1)*N
        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                            // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;              // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (s=0,nb1i=0ull,k=0;k<8;k++) {                                      // packs non-zero nb indices in first up to 8*4 bits
            gols=gol[nb[k]];                                                   // whether neighbor is alive
            s += gols;                                                         // s is number of live nbs
            nb1i = (nb1i << (gols<<2)) + (gols*k);                             // nb1i is packed list of live neighbour indices
        }
        smid = s>1 && s<7; s2 = s-2;                                           // s in mid-range for possible lut rule
        statflag = newgene = survive = birth = 0ull;
        if (smid) {
            s2or3 = (s>>2) ? 0ull : (s>>1);                                    // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
            gols = s2or3 ? (gol[ij] ? 1ull : (s&1ull ? 1ull : 0ull )) : 0ull;  // GoL calculation next state for non-genetic gol plane
            rulemodij = (rulemod&0x2) ? (ij>=(N2>>1) ? 1 : 0) : (rulemod&0x1); // if rulemod bit 1 is on then split into half planes with/without mod
            if (rulemodij) {
                overwrite = s ? (overwritemask>>(s-1))&0x1ull : 0ull;          // allow birth to overwrite occupied cell = survival in GoL
                overwrite = overwrite | (~gol[ij] & 0x1ull);                   // either central cell is empty or overwrite bit set is required for birth
                if (gol[ij]) survive = (smask>>sumoffs[s2])&summasks[s2] ? 1ull : 0ull;
                if (overwrite) birth = (bmask>>sumoffs[s2])&summasks[s2] ? 1ull : 0ull;
                /* if(s==3) fprintf(stderr,"lut_canon step %d ij %d s %d gol[ij] %llx overwrite %llx survive %llx birth %llx\n",
                                                     totsteps,ij,s,gol[ij],overwrite,survive,birth); */
                if (survive|birth) {
                    for (k=0,nbmask=0;k<8;k++) nbmask |= (gol[nb[k]]<<k);      // constuct mask of live bits for 8 neighbours
                    kch = selectdifft(s, nbmask, &crot, &kodd, &nsame);        // find the canonical rotation index of the live neighbour configuration
                    survive= (smask>>(sumoffs[s2]+crot))&0x1ull;               // refine decisions for specific canonical rotation configuration
                    birth  = (bmask>>(sumoffs[s2]+crot))&0x1ull;               // only allowed if birth/survivemask permits (ask this before consulting genes)
                    if (survive|birth) {                                       // complete determination of birth or survival
                        if (selection==13) {                                   // modular gene encoding one of up to two 5-bit subsets of crot (e.g. for s=3,4 5+2 and 5+5 resp.)
                            survive=birth=0ull;
                            crot5 = crot > 4 ? 1 : 0;                          // crot is in range 0-9
                            crotmod5 = crot-5*crot5;                           // the 5 bits indexed by crotmod5 are stored in two places 0-2 in 8-bit word and 3-4 in pairs bits 48+
                            pat = (s<<4)|(crot5<<3)|(0x1<<((crotmod5>2)?crotmod5-3:crotmod5));
                            genecode = 0ull;
                            for (k=0;k<8;k++) {                                // decodes genes with variable length encoding only for current s,crot
                                if (gol[nb[k]]) {                              // combine information from genes of all live neighbours
                                    if ((~survivalgene & gol[ij]) | overwrite) {
                                        genecode = golg[nb[k]];
                                        if(crotmod5 > 2) for (genecode1=0ull,k1=0;k1<6;k1++) genecode1 |= ( (genecode>>(48+(k1<<1))) & (0x1ull<<(crotmod5-3)) ) << (k1<<3);
                                        else genecode1 = genecode & (0x010101010101 << crotmod5);     // prepare lookup in all 6 modules on gene : in this case crot subset bits in right place
                                        genecode &= 0xf8f8f8f8f8f8; //
                                        genecode |= genecode1;                 // replace lowest 3 bits of 8-bit part of 10-bit modules with appropriate crot subset bits
                                    }
                                    if (~survivalgene & gol[ij]) {
                                        PATTERN8(genecode, pat, found);
                                        survive |= found? 1ull : 0ull;         //final decision for survival?
                                    }
                                    if (overwrite) {
                                        PATTERN8(genecode, (0x80|pat), found);
                                        birth |= found? 1ull : 0ull;           //final decision for birth?
                                    }
                                }
                            }
                            /* fprintf(stderr,"lut_canon step %d ij %d s %d gol[ij] %llx crot %d (%d,%d) pat %x genecode %llx overwrite %llx survive %llx birth %llx\n",
                                                     totsteps,ij,s,gol[ij],crot,crot5,crotmod5,pat,genecode,overwrite,survive,birth); */
                            // if (overwrite && !gol[ij] && s==3) fprintf(stderr,"overwrite and golij!!!!!!!!!!!!!\n");
                            if (survivalgene & gol[ij]) {                     // survival determined by central gene in this case
                                genecode = golg[ij];
                                if(crotmod5 > 2) for (genecode1=0ull,k1=0;k1<6;k1++) genecode1 |= ((genecode>>(48+(k1<<1)))&(0x1<<(crotmod5-3)))<<(k1<<3);
                                else genecode1 = genecode & (0x010101010101 << crotmod5);
                                genecode &= 0xf8f8f8f8f8f8;
                                genecode |= genecode1;
                                PATTERN8(genecode, pat, found); //final decision for survival?
                                survive = found? 1ull : 0ull;
                            }
                        }
                        else {                                                 // selection == 12
                            if(repscheme&R_1_nb_OR_AND)
                                for (genecode=0ull,k=0;k<8;k++)                     // decodes genes with fixed length encoding by OR
                                    genecode |= (gol[nb[k]]?golg[nb[k]]:0ull);      // OR of live neighbours encodes birth rule & survival rule
                            else
                                for (genecode=~0ull,k=0;k<8;k++)                    // decodes genes with fixed length encoding by AND
                                    genecode &= (gol[nb[k]]?golg[nb[k]]:~0ull);     // AND of live neighbours encodes birth rule & survival rule
                            if(survivalgene) genecode = (genecode&(0xffffffffull<<32)) | (golg[ij]&0xffffffffull);  // if central gene determines survival
                            genecode&=(bmask<<32)|smask;                       // actually no longer needed since test done above
                            if (gol[ij]) survive |= (genecode>>(sumoffs[s2]+crot))&0x1ull;
                            if (overwrite) birth &= (genecode>>(32+sumoffs[s2]+crot))&0x1ull;
                        }
                    }
                }
            }
            else {
                survive = s2or3;
                birth = s2or3&s&0x1&~gol[ij];
            }
        }
        else survive = birth = s2or3 = gols = 0ull;

        if(birth) {
            if (repscheme & R_7_random_resln) {
                RAND128P(randnr);                                               // inline exp so compiler recognizes auto-vec,
                kch = ((randnr>>32)&0xffff) % s;                                // choose random in this option only
                newgene = golg[nb[(nb1i>>(kch<<2))&0x7]];
            }
            else {
                nbest=selectone_of_s(s,nb1i,nb,golg,&birth,&newgene,&nbmask);   // selection scheme depends on repscheme parameter
                if(nbest>1) {
                    kch=selectdifft(nbest,nbmask, &crot, &kodd, &nsame);
                    RAND128P(randnr);                                               // inline exp so compiler recognizes auto-vec,
                    if(nsame) newgene = disambiguate(kch, nb1i, nb, golg, nsame, &birth, randnr); // restore symmetry via one of 8 repscheme options
                    else newgene = golg[nb[kch]];
                }
            }
            if (birth) {                                                        // renewed query as possibly updated in disambiguate
                statflag |= F_birth;
                r2=1ull;                                                        // compute random events for single bit mutation, as well as mutation position nmut
                ancestor = newgene;
                while (r2) {                                                    // allow multiple point mutations
                    RAND128P(randnr);                                           // inline exp so compiler recognizes auto-vec,
                    r2 = randprob(pmutmask,(unsigned int) randnr);
                    nmut = (randnr >> 56) & 0x3f;                               // choose mutation position for length 64 gene : from bits 56:61 of randnr
                    newgene = newgene ^ (r2<<nmut);                             // introduce single mutation with probability pmut = probmut
                    statflag = statflag | F_mutation;
                }
                newgol[ij]  = 1ull;
                newgolg[ij] =  newgene;                                         // if birth then newgene
                if(gol[ij]) hashdeletegene(golg[ij],"step %d hash delete error 1 in update_lut_sum, gene %llx not stored\n");
                hashaddgene(newgene,ancestor);
            }
            else {                                                              // no birth, stays empty
                newgol[ij] = gol[ij];
                newgolg[ij] = golg[ij];
            }
        }
        else if(gol[ij]) {                                                      // death/survival   (careful! genecode changed during this processing)
            if(survive) {
                statflag |= F_survival;
                newgol[ij]  = gol[ij];                                          // new game of life cell value same as old
                newgolg[ij] = golg[ij];                                         // gene stays same
            }
            else {                                                              // gene bit = 0 => death
                statflag |= F_death;
                newgol[ij]  = 0ull;                                             // new game of life cell value dead
                newgolg[ij] = 0ull;                                             // gene dies
                hashdeletegene(golg[ij],"step %d hash delete error 2 in update, gene %llx not stored\n");
            }
        }
        else{                                                                   // no birth on empty cell
            newgol[ij] = gol[ij];
            newgolg[ij] = golg[ij];
        }
        
        if(newgol[ij]!=gols) {
            statflag |= F_notgolrul;
            if(newgol[ij]) statflag |= F_nongolchg;
        }

        if(gol[ij]) statflag |= F_golstate;                                 // this is the last gol state, not the new state
        newgolgstats[ij] = statflag;
    }  // end for ij

    if(randomsoup) random_soup(gol,golg,newgol,newgolg);
    if(vscrolling) v_scroll(newgol,newgolg);
    if (colorfunction==8) packandcompare(newgol,working,golmix);

    for (ij=0; ij<N2; ij++) {       // complete missing hash table records of extinction and activities
        if(gol[ij]) hashgeneextinction(golg[ij],"hash extinction storage error in update, gene %llx not stored\n");
        if(newgol[ij]) hashgeneactivity(newgolg[ij],"hash activity storage error in update, gene %llx not stored\n");
    }
    if(quadtree) qimage = quadimage(newgol,N,0);
}
//------------------------------------------------------- update_lut_2Dsym -----------------------------------------------------------------------------------
void update_lut_2D_sym(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]){     // selection models 14,15
/*
    all different configurations under the standard 2D 4-rotation and 4-reflection symmetries are distinguished
    i.e. by number of ones and edge-corner differences and additional distinctions in arrangement
    for the different s values                     0  1  2  3  4  5  6  7  8
    there are a nr of distinct configurations      1  2  6 10 13 10  6  2  1   total 51
    if we exclude the cases s=0,1,7,8 which can create growth artefacts, then there are still 45 bits
    if we focus on cases s=3,4,5 there are still one too many, i.e. 33 cases, for investigation in a single integer genome
    instead, we investigate here the most interesting lower s domain s = 0,1,2,3,4 in full detail
    [     alternative approaches would be to:
        (i) investigate separately the configurations s=4-8 which requires 32 bits
        (ii) investigate s=2,3,5,6, excluding 4, with 32 bits for birth and survival
        (iii) code up the s=3,4,5 minus one : with the software specifying one omitted configuration
        (iv) include undifferentiated s=5 instead of 0 : 32 bits
        (v) include undifferentiated s=1,5,6 instead of 0,1  ie s=1-6 : 32 bits
        (vi) split the s-range allowed for birth and survival: S s=1-4 B s=3-5 ie 31 bits for survival and 33 for birth]
 
                                                                                                                */
    int s, se, slow, s2or3, k, k1, kodd, coff, nmut, crot, found, coffdiv6, coffmod6, pat, nbest, nsame;
    uint64_t survive,birth,overwrite,survivalgene,smask,bmask,rulemodij;
    static uint64_t summasks[5] = {0x1ull,0x3ull,0x3full,0x3ffull,0x1fffull};  // masks for s= 0,1,2,3,4 with nr cases 1,2,6,10,13
    static int sumoffs[9] = {0,1,3,9,19,32,42,48,50};                                      // cumulative offsets to start of coding region for s = 0,1,2,3,4,5,6,7,8
    static int csumoffs[9] = {0,2,4,12,26,46,60,68,2};                                     // start of indexing in confoffs for crot,kodd lookup for s = 0,1,2,3,4,5,6,7,8
    static unsigned char confoffs[2+2+8+14+20+14+8+2+2] = {0,0, 0,0, 0,0,1,2,3,3,4,5, 0,1,2,3,3,2,4,5,6,7,4,5,8,9, 0,0,1,2,3,4,1,2,5,6,7,8,9,9,10,10,8,7,11,12,
                                        0,1,2,3,3,2,4,5,6,7,4,5,8,9, 0,0,1,2,3,3,4,5, 0,0, 0,0}; // look up for crot*2+kodd to gene bit offset
    static unsigned char lut[256];
    int nb[8], ij, i, j, jp1, jm1, ip1, im1;
    unsigned int kch=0;
    uint64_t genecode, genecode1, gols, nb1i, nbmask, randnr, r2;
    uint64_t newgene, ancestor;
    uint64_t statflag;
    
    static int first = 1;
    
    if (first) {
        first = 0;
        for (nbmask=0;nbmask<256;nbmask++) {
            POPCOUNT64C(nbmask,s);
            kch = selectdifft(s, nbmask, &crot, &kodd, &nsame);
            coff = confoffs[csumoffs[s]+(crot<<1)+kodd];
            lut[nbmask]=coff;
            fprintf(stderr,"LUT nbmask %2llx s %1d goff %2d soff %2d coff %2d crot %2d kodd %1d nsame %2d\n",nbmask,s,sumoffs[s]+coff,sumoffs[s],coff,crot,kodd,nsame);
        }
    }
    
    canonical = repscheme & R_2_canonical_nb;                                  // set global choice of canonical rotation bit choice for selectdifftx
    survivalgene = repscheme & R_0_survivalgene;                               // gene determining survival is 1: central gene 2: determined by neighbours
    smask = (uint64_t) survivalmask;                                           // 32 bits of survivalmask used to limit space of rules, convert to 64 bit masks for efficient usage here
    bmask = (uint64_t) birthmask;                                              // 32 bits of birthmask

    for (ij=0; ij<N2; ij++) {                                                  // loop over all sites of 2D torus with side length N
        i = ij & Nmask;  j = ij >> log2N;                                      // row & column
        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                      // toroidal (j+1)*N and (j-1)*N
        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                            // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;              // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (s=se=0,nb1i=0ull,k=0;k<8;k++) {                                   // packs non-zero nb indices in first up to 8*4 bits
            gols=gol[nb[k]];                                                   // whether neighbor is alive
            s += gols;                                                         // s is number of live nbs
            se += k&0x1&gols;                                                  // se is number of edge-centred live neighbours (odd k)
            nb1i = (nb1i << (gols<<2)) + (gols*k);                             // nb1i is packed list of live neighbour indices
        }
        slow = s<5;                                                            // s in low-range 0-4 for possible lut rule
        statflag = newgene = survive = birth = 0ull;
        if (slow) {
            s2or3 = (s>>2) ? 0ull : (s>>1);                                    // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
            gols = s2or3 ? (gol[ij] ? 1ull : (s&1ull ? 1ull : 0ull )) : 0ull;  // GoL calculation next state for non-genetic gol plane
            rulemodij = (rulemod&0x2) ? (ij>=(N2>>1) ? 1 : 0) : (rulemod&0x1); // if rulemod bit 1 is on then split into half planes with/without mod
            if (rulemodij) {
                overwrite = s ? overwritemask&(0x1ull<<(s-1)) : 0;
                overwrite = (overwrite || !gol[ij]) ? 1ull : 0ull;
                if (gol[ij]) survive = (smask>>sumoffs[s])&summasks[s] ? 1ull : 0ull;
                if (overwrite) birth = (bmask>>sumoffs[s])&summasks[s] ? 1ull : 0ull;
                if (survive|birth) {
                    for (k=0,nbmask=0;k<8;k++) nbmask |= (gol[nb[k]]<<k);      // constuct mask of live bits for 8 neighbours
                    kch = selectdifft(s, nbmask, &crot, &kodd, &nsame);        // find the canonical rotation index and odd/even offset of the live neighbour configuration
                    coff = confoffs[csumoffs[s]+(crot<<1)+kodd];
                    coffdiv6 = coff >= 12 ? 2 : (coff >=6 ? 1 : 0);            // coff is in range 0-12 coffdiv6 0,1,2
                    coffmod6 = coff-6*coffdiv6;                                // the 6 bits indexed by coffmod6 are stored in two places 0-1 in 8-bit word and 2-5 in quartet bits 40+
                    pat = (s<<4)|(coffdiv6<<2)|0x1;
                    genecode = 0ull;
                    survive= (smask>>(sumoffs[s]+coff))&0x1ull;                           // refine decisions for specific canonical rotation configuration
                    birth  = (bmask>>(sumoffs[s]+coff))&0x1ull;                           // only allowed if birth/survivemask permits (ask this before consulting genes)
                    if (survive|birth) {                                       // complete determination of birth or survival
                        if (selection==15) {                                   // selection == 15 variable length encoding
                            survive=birth=0ull;
                            for (k=0;k<8;k++) {                                // decodes genes with variable length encoding only for current s,crot
                                if (gol[nb[k]]) {                              // combine information from genes of all live neighbours
                                    if ((!survivalgene && gol[ij]) || overwrite) {
                                        genecode = golg[nb[k]];
                                        if(coffmod6 > 1) for (genecode1=0ull,k1=0;k1<5;k1++) genecode1 |= ( (genecode>>(40+(k1<<2)+coffmod6-2)) & 0x1ull ) << (k1<<3);
                                        else genecode1 = (genecode >> coffmod6) & 0x0101010101; // prepare lookup in all 5 modules on gene : in this case coff subset bits in right place
                                        genecode &= 0xfcfcfcfcfc; //
                                        genecode |= genecode1;                 // replace lowest 2 bits of 8-bit part of 12-bit modules with appropriate coff subset bits
                                    }
                                    if (!survivalgene && gol[ij]) {
                                        PATTERN8(genecode, pat, found); //final decision for survival?  NB all 0 seq encodes survival with 0 live nbs
                                        survive |= found? 1ull : 0ull;
                                    }
                                    if (overwrite) {
                                        PATTERN8(genecode, (0x80|pat), found); //final decision for birth?
                                        birth |= found? 1ull : 0ull;
                                    }
                                }
                            }
                            if (survivalgene && gol[ij]) {                     // survival determined by central gene in this case
                                genecode = golg[ij];
                                if(coffmod6 > 1) for (genecode1=0ull,k1=0;k1<5;k1++) genecode1 |= ( (genecode>>(40+(k1<<2)+coffmod6-2)) & 0x1ull ) << (k1<<3);
                                else genecode1 = (genecode >> coffmod6) & 0x0101010101; // prepare lookup in all 5 modules on gene : in this case coff subset bits in right place
                                genecode &= 0xfcfcfcfcfc; //
                                genecode |= genecode1;                 // replace lowest 2 bits of 8-bit part of 12-bit modules with appropriate coff subset bits
                                PATTERN8(genecode, pat, found);      //final decision for survival?  NB all 0 seq encodes survival with 0 live nbs
                                survive |= found? 1ull : 0ull;
                            }
                        }
                        else {                                                 // selection == 14
                            if(repscheme&R_1_nb_OR_AND)
                                for (genecode=0ull,k=0;k<8;k++)                     // decodes genes with fixed length encoding by OR
                                    genecode |= (gol[nb[k]]?golg[nb[k]]:0ull);      // OR of live neighbours encodes birth rule & survival rule
                            else
                                for (genecode=~0ull,k=0;k<8;k++)                    // decodes genes with fixed length encoding by AND
                                    genecode &= (gol[nb[k]]?golg[nb[k]]:~0ull);     // AND of live neighbours encodes birth rule & survival rule
                            if(survivalgene) genecode = (genecode&(0xffffffffull<<32)) | (golg[ij]&0xffffffffull);  // if central gene determines survival
                            genecode&=(bmask<<32)|smask;                       // actually no longer needed since test done above
                            if (gol[ij]) survive |= (genecode>>(sumoffs[s]+coff))&0x1ull;
                            if (overwrite) birth &= (genecode>>(32+sumoffs[s]+coff))&0x1ull;
                        }
                    }
                }
            }
            else {
                survive = s2or3;
                birth = s2or3&s&0x1&~gol[ij];
            }
        }
        else survive = birth = s2or3 = gols = 0ull;

        if(birth) {
            if (repscheme & R_7_random_resln) {
                RAND128P(randnr);                                               // inline exp so compiler recognizes auto-vec,
                kch = ((randnr>>32)&0xffff) % s;                            // choose random in this option only
                newgene = golg[nb[(nb1i>>(kch<<2))&0x7]];
            }
            else {
                nbest=selectone_of_s(s,nb1i,nb,golg,&birth,&newgene,&nbmask); // selection scheme depends on repscheme parameter
                if(nbest>1) {
                    kch=selectdifft(nbest,nbmask, &crot, &kodd, &nsame);
                    RAND128P(randnr);                                               // inline exp so compiler recognizes auto-vec,
                    if(nsame) newgene = disambiguate(kch, nb1i, nb, golg, nsame, &birth, randnr); // restore symmetry via one of 8 repscheme options
                    else newgene = golg[nb[kch]];
                }
            }
            if (birth) {                                                    // renewed query as possibly updated in disambiguate
                statflag |= F_birth;
                r2=1ull;                                                        // compute random events for single bit mutation, as well as mutation position nmut
                ancestor = newgene;
                while (r2) {                                                    // allow multiple point mutations
                    RAND128P(randnr);                                           // inline exp so compiler recognizes auto-vec,
                    r2 = randprob(pmutmask,(unsigned int) randnr);
                    nmut = (randnr >> 56) & 0x3f;                               // choose mutation position for length 64 gene : from bits 56:61 of randnr
                    newgene = newgene ^ (r2<<nmut);                             // introduce single mutation with probability pmut = probmut
                    statflag = statflag | F_mutation;
                }
                newgol[ij]  = 1ull;
                newgolg[ij] =  newgene;                                     // if birth then newgene
                if(gol[ij]) hashdeletegene(golg[ij],"step %d hash delete error 1 in update_lut_sum, gene %llx not stored\n");
                    hashaddgene(newgene,ancestor);
            }
            else {                                                          // no birth, stays empty
                newgol[ij] = gol[ij];
                newgolg[ij] = golg[ij];
            }
        }
        else if(gol[ij]) {                                                      // death/survival   (careful! genecode changed during this processing)
            if(survive) {
                statflag |= F_survival;
                newgol[ij]  = gol[ij];                                          // new game of life cell value same as old
                newgolg[ij] = golg[ij];                                         // gene stays same
            }
            else {                                                              // gene bit = 0 => death
                statflag |= F_death;
                newgol[ij]  = 0ull;                                             // new game of life cell value dead
                newgolg[ij] = 0ull;                                             // gene dies
                hashdeletegene(golg[ij],"step %d hash delete error 2 in update, gene %llx not stored\n");
            }
        }
        else{                                                                   // no birth on empty cell
            newgol[ij] = gol[ij];
            newgolg[ij] = golg[ij];
        }
        
        if(newgol[ij]!=gols) {
            statflag |= F_notgolrul;
            if(newgol[ij]) statflag |= F_nongolchg;
        }

        if(gol[ij]) statflag |= F_golstate;                                 // this is the last gol state, not the new state
        newgolgstats[ij] = statflag;
    }  // end for ij

    if(randomsoup) random_soup(gol,golg,newgol,newgolg);
    if(vscrolling) v_scroll(newgol,newgolg);
    if (colorfunction==8) packandcompare(newgol,working,golmix);

    for (ij=0; ij<N2; ij++) {       // complete missing hash table records of extinction and activities
        if(gol[ij]) hashgeneextinction(golg[ij],"hash extinction storage error in update, gene %llx not stored\n");
        if(newgol[ij]) hashgeneactivity(newgolg[ij],"hash activity storage error in update, gene %llx not stored\n");
    }
    if(quadtree) qimage = quadimage(newgol,N,0);
}
//------------------------------------------------------- update_gol16 -----------------------------------------------------------------------------------
void update_gol16(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]) {     // selection models 16-19
// genes specify on which of 16 planes the gol states at that site are visible: specified by 16 bits 0,4,8,... in gene
// when multiple copy processes are active for gol states, the one with the lowest bit position is applied to most difft gene copying
// gol states looked up for neighbors are the OR of the current plane gol value with gol values on planes specified by the gene for repscheme&1 on
// since extra birth is achievewd by coupling, we can optionally (repscheme&2) not allow survival to balance growth
// option 1 (sel 10) genes also specify a unique copy plane for a gene: it can only be copied by a rule active on this plane
// option 2 (sel 11) coupling is to nearby planes (2 or 4 depending on repscheme&4) with bitmask specified by genes
// alternatively (repscheme&8) copy plane could be lowest non zero bit in gene at pos m for which m%4 == 1 (if none then gene is not copied)
// repscheme control bits .........................................................................................
// RP_0_or_current   0x1        1: do or between coupling plane and current plane(s) 0: just take coupling plane(s)
// RP_1_c4planes     0x2        1: 4 nearest planes coupled 0: 2 nearest planes coupled
// RP_2_survive3     0x4        1: survival for 3 live-nbs 0: no survival for 3 live-nbs
// RP_3_survive2     0x8        1: survival for 2 live-nbs 0: no survival for 2 live-nbs
// RP_4_grow_planes  0x10       1: grow from plane 0 0: chose coupling to any plane

    int crot,kodd,c4planes;
    unsigned int ij,ij1,k,kch,p,p1,nmut,debcnt,mask,np,npmask;
    uint64_t s,sm,sm3,s3sm3,su,s3,sg3,s2or3,nbmask,ancestor,newgene,golsh,pmask;
    uint64_t *golm;
    //uint64_t s0,sm;
    uint64_t randnr, rand2, statflag;
    int nb1x[8] = {-1,0,1,1,1,0,-1,-1};
    int nb1y[8] = {-1,-1,-1,0,1,1,1,0};
    int npmasks[17] = {1,1,3,3,7,7,7,7,0xf,0xf,0xf,0xf,0xf,0xf,0xf,0xf,0xf};
    const uint64_t r1= 0x1111111111111111;
    const uint64_t rc= 0xcccccccccccccccc;
    
    np = NbP;                                                               // number of planes used, 16 if 0
    npmask = npmasks[np-1];                                                 // mask for other plane index
    pmask = np==16 ? 0xffffffffffffffff : (1ull<<(np<<2))-1;                // mask for bits included in np planes
    
    if(!(totsteps%10)) fprintf(stderr,"iteration step %d\r",totsteps);
    if(rulemod) {
      if (selection >= 16 && selection <=17) {                               // gene codes for one secondary plane for each plane  (OLD 10)
        for (ij=0; ij<N2; ij++) {                                            // loop over all sites of 2D torus with side length N
            for (p=0,golsh=0ull;p<np;p++) {                                  // 16 different planes, one every 4 bits in gol,golg
                p1 = (golg[ij]>>(p<<2))&npmask;
                p1 &= (repscheme&RP_4_grow_planes) ? npmasks[p+1] : 0xf;     // 1st implementation of Norman's idea to grow plane community
                golsh |= ((gol[ij]>>(p1<<2))&1ull)<<(p<<2);                  // collect active bits from 2nd planes pointed to by gene
            }
            golsh &= pmask;
            if (repscheme&RP_0_or_current)   golmix[ij] = golsh|gol[ij];     // bits active from secondary plane or primary plane
            else golmix[ij] = golsh;                                         // bits active from secondary plane only
        }
      }
      else { // selection == 18,19                                           // gene codes for 4(2) nearest plane neighbour masks
        if (repscheme&RP_1_c4planes) c4planes = 1;                           // decision 4 or 2 based on repscheme bit
        else c4planes=0;
        for (ij=0; ij<N2; ij++) {                                            // loop over all sites of 2D torus with side length N
            for (p=0,golsh=0ull;p<np;p++) {                                  // 16 different planes, one every 4 bits in gol,golg
                mask = (golg[ij]>>(p<<2))&0xf;                               // bit mask indicating which nearby planes visible from plane
                for (k=1-c4planes;k<3+c4planes;k++) {                        // 2 or 4 nearest planes depending on repscheme bit 2
                    if((mask>>k)&0x1) {
                        p1 = (p+k+(k>>1)-2)&npmask;
                        golsh ^= ((gol[ij]>>(p1<<2))&0x1ull)<<(p<<2);        // nearby planes p1=p+{-2,-1,1,2} for k 0,1,2,3   or p+{-1,1} for k 1,2
                    }
                }
            }
            golsh &= pmask;
            if (repscheme&RP_0_or_current)   golmix[ij] = golsh|gol[ij];     // bits active from secondary plane or primary plane
            else golmix[ij] = golsh;                                         // bits active from secondary plane only
        }
      }
    }
    else { // no gol rule modification, independent planes
         for (ij=0; ij<N2; ij++) golmix[ij] = gol[ij];                       // loop over all sites of 2D torus with side length N
    }
    // for (ij=0; ij<N2; ij++) if ((golmix[ij]&r1) != golmix[ij]) fprintf(stderr,"error in golmix calculation %llx",golmix[ij]);
    
    for (ij=0; ij<N2; ij++) {                                                // loop over all sites of 2D torus with side length N
        for(k=0,s=0ull,sm=0ull;k<8;k++) {                                    // compute no of live neighbours using golmix (coupled plane neighbours)
            ij1=deltaxy(ij,nb1x[k],nb1y[k]);
            sm +=golmix[ij1];
            s +=gol[ij1];
        }
        sm3 = (~sm>>3)&(~sm>>2)&(sm>>1)&sm&r1;                               // sum of golmix bits is exactly 3, calculated for each plane
        su = s&rc;                                                           // upper 2 bits (3,2) of sum of live neighbours : non zero if sum 4-8
        sg3 = (((su>>1)|su)>>2) & r1;                                        // for each plane 1 if sum is gt 3
        s2or3 = (~sg3)&(s>>1)&r1;                                            // sum is 2 or 3 for each plane, ie not gt 3 and bit 2 is 1
        s3 = s2or3&s;                                                        // sum is 3 (s2or3 and bit 0 is 1) calculated for each plane
        s3sm3 = s3|sm3;
        newgol[ij]=s2or3&s3sm3;                                              // only birth so far
        if (repscheme&RP_2_survive3)                                         // add survival depending on repscheme
            newgol[ij]|=gol[ij]&s3;
        if (repscheme&RP_3_survive2)
            newgol[ij]|=gol[ij]&(s2or3&~s3);
        statflag = 0ull;
        //if((golmix[ij]^gol[ij])&s2or3) statflag |= F_notgolrul;            // not gol rule in at least one plane
        if((s3^sm3)&s2or3) statflag |= F_notgolrul;                          // not gol rule in at least one plane
        if (s3sm3) {                                                         // birth in at least one plane, implement for lowest such birth plane
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
                kch=selectdifft3(nbmask, &crot, &kodd);
                ij1=deltaxy(ij,nb1x[kch],nb1y[kch]);
                newgene = golg[ij1];
                RAND128P(randnr);                                            // inline exp so compiler recognizes auto-vec,
                rand2 = randprob(pmutmask, (unsigned int) randnr);
                nmut = (randnr >> 56) & (0x3|(npmask<<2));                   // choose mutation pos for length np*4 gene : from bits 56:61 of randnr
                ancestor = newgene;
                newgene = newgene ^ (rand2<<nmut);                           // introduce single mutation with probability pmut = probmut
                newgolg[ij]=newgene;
                statflag = statflag | F_birth;
                if (rand2) statflag = statflag | F_mutation;
                newgolgstats[ij] = statflag;
                hashreplacegene(golg[ij],newgene,ancestor,"step %d hash storage error 1 in update_gol16, gene %llx not stored\n");
                p=np;                                                        // break from loop at first birth
            } //if 3 live nbs at p
        } // for p
      }  // if (s3)
      else {
        newgolg[ij]=gene0;
        if(golg[ij]!= gene0) {
            newgene= gene0;                                                  // default is relax to uncoupled gene
            hashreplacegene(golg[ij],newgene,rootgene,"step %d hash storage error 1 in update_gol16, gene %llx not stored\n");
        }
      } // end else (s3)
    } // for ij
    
    if(randomsoup) random_soup(gol,golg,newgol,newgolg);
    if(vscrolling) v_scroll(newgol,newgolg);
    for (ij=0; ij<N2; ij++) {  // complete missing hash table records including activities, NB all sites have genes in gol16 but don't record activity for zerogene
        hashgeneextinction(golg[ij],"hash storage error 4 in update_gol16, gene %llx not stored\n");
        if(newgolg[ij]!=gene0) hashgeneactivity(newgolg[ij],"hash storage error 5 in update_gol16, gene %llx not stored\n");
    }
}
//------------------------------------------------------- update_gol64 -----------------------------------------------------------------------------------
void update_gol64(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]) {  // selection model 20+, routine for 64x-gol packed update with 1 gene plane
    unsigned int ij,ij1,k,k1,kch,kb,sum,nmut,nbmask,birthinitplane,plane;
    uint64_t golij,gs,newgs,newgs3d,sums00,sums10,sums01,sums11,sums02,sums12,sums16[4],sums3D[4];
    uint64_t s,s3,s4,s5,s6,s7,s567,sgol,sm,sh,sge8,plcodingmask,golgij,newgene,ancestor,b,b3,b3d,genemask;
    // uint64_t s2shR,s2shL,g;
    uint64_t randnr, rand2, statflag;
    int crot,kodd;
    int nb1x[9] = {-1,0,1,1,1,0,-1,-1,0};
    int nb1y[9] = {-1,-1,-1,0,1,1,1,0,0};
    const uint64_t r1= 0x1111111111111111;
    const uint64_t r3= 0x3333333333333333;
    const uint64_t r5= 0x5555555555555555;
    const uint64_t r7= 0x7777777777777777;
    const uint64_t r8= 0x8888888888888888;
    const uint64_t ra= 0xaaaaaaaaaaaaaaaa;
    const uint64_t rc= 0xcccccccccccccccc;
    
    canonical = repscheme & R_2_canonical_nb;
    plcodingmask= (NbP==64) ? ~0ull : (1ull<<NbP)-1;                            // note that for number of planes less than 64, the 3D calculation is not periodic closed. PROBLEM! *******
    if(!(totsteps%10)) fprintf(stderr,"iteration step %d\r",totsteps);
    
    for (ij=0; ij<N2; ij++) {                                                   // loop over all sites of 2D torus with side length N
        golij = gol[ij];
        if(golij && golg[ij]==0ull) fprintf(stderr,"First step %d ij %d golij non zero %llx but golgij zero\n",totsteps,ij,golij);
        sums00 = sums10 = sums01 = sums11 = sums02 = sums12 = 0ull;             // assemble 4 bit sums of planar 9-neighborhoods
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
        for(k=6;k<9;k++) {
            gs=gol[deltaxy(ij,nb1x[k],nb1y[k])];
            sums02 += gs&r5;
            sums12 += (gs&ra)>>1;
        }
        sums16[0] = (sums00&r3)+(sums01&r3)+(sums02&r3);
        sums16[1] = (sums10&r3)+(sums11&r3)+(sums12&r3);
        sums16[2] = ((sums00&rc)>>2)+((sums01&rc)>>2)+((sums02&rc)>>2);
        sums16[3] = ((sums10&rc)>>2)+((sums11&rc)>>2)+((sums12&rc)>>2);
        newgs=0ull;b3=b=0ull;
        for (k=0;k<4;k++) {                                                     // 2D GoL calculation "egocentrically" see Wikipedia GoL page for term
            s  = sums16[k];
            s4 =(~s>>3)&((s>>1)^~s);                                            // do common part of s3,s4 calculation first: bit3 0 and bit0=bit1
            s3 = s4 & (~s>>2)&s&r1;
            s4 = s4 & (s>>2)&~s&r1;
            sgol=s3|(s4&(golij>>k)&r1);
            newgs|=sgol<<k;
            // s2 |= (s3&(golij>>k)&r1)<<k;
            b3 |= (s3&(~golij>>k)&r1)<<k;
        }
        b= newgs&~golij;
        // if (b!=b3) fprintf(stderr,"in update_gol64 something went wrong with 2D-GoL birth process b=%llx != b3=%llx\n",b,b3);
        if (rulemod) {
            for (k=0;k<4;k++) {
                sm = s = sums16[k]-((golij>>k)&r1);                             // the sum for each central plane should not include the central site
                sge8 = s&r8;                                                    // bit 3 is 1 if the corresponding plane sum is >= 8
                sm &= r7;
                k1=(k+1)&0x3;
                sh = sums16[k1];                                                // now add in sum of all 9 sites from next plane above
                sh = k1 ? sh : (sh>>4)|(sh<<60);
                sm += sh&r7;
                sge8 |= (sm|sh)&r8;
                sm &= r7;
                k1=(k-1)&0x3;
                sh = sums16[k1];                                                // now add in sum of all 9 sites from next plane below
                sh = k ? sh : (sh<<4)|(sh>>60);
                sm += sh&r7;
                sums3D[k] = sm|sge8|(sh&r8);                                    // sums3D are correct modulo 8 for 26 nbs and the 4th bit determines if sum is >=8.
            }
            // determine if values are 567 for survival or 6 for birth (excluding central cell): Carter Bay's 3D life rule 5766.
            for (s5=s6=s7=s567=0ull,k=0;k<4;k++) {
                s  = sums3D[k];                                                 // sums in 3D calculated above
                sm = (~s>>3)&(s>>2);                                            // partial calculation of common bits between s6 and s567
                sh = sm&(s>>1)&~s&r1;                                           // contribution to s6 from k
                s6 |= sh<<k;                                                    // assembled bits for all 64 planes of whether sum=6
                s567 |= (sh|(sm&s&r1))<<k;                                      // assembled bits for all 64 planes of whether sum=5,6 or 7
                s5 |= (sm&~(s>>1)&s&r1)<<k;                                     // needed for repscheme coupling variants not for standard 3D
                s7 |= (sm&(s>>1)&s&r1)<<k;                                      // needed for repscheme coupling variants not for standard 3D
            }
            b3d = s6 & ~golij;
            if      (repscheme&0x10) newgs3d = s6|(s5&golij);                   // new state is modified 3D GoL next state with survival only for 5 or 6 not 7
            else if (repscheme&0x20) newgs3d = s6|(s7&golij);                   // new state is modified 3D GoL next state with survival only for 6 or 7 not 5
            else if (repscheme&0x40) newgs3d = s6;                              // new state is modified 3D GoL next state with survival only for 6 not 7 not 5
            else                     newgs3d = s6|(s567&golij);                 // new state is 3D GoL next state
        }
        else b3d=newgs3d=0ull;

        golgij=golg[ij];
        if(rulemod) {
            if(repscheme&0x1) {                                                 // turn coupling on: 2D plus 3D when gene specifies the plane
                if(repscheme&0x2) {                                             // with this repscheme parameter, only use birth not survival from 2D GoL
                    newgs = b3;
                }
                //genemask = 1ull<<(golgij&0x3full);                            // gene allows coupling in plane specified by last 6 bits of gene
                genemask = golgij;                                              // gene allows coupling in planes where a 1 is in genome only
                b3d&=genemask;
                newgs3d&=genemask;
                b = b3|b3d;
                newgol[ij]=(newgs|newgs3d)&plcodingmask;
            }
            else {                                                              // genes have no effect on pure 3D GoL for repscheme = 0
                b = b3d;
                b3=0ull;
                newgol[ij]=newgs3d&plcodingmask;
            }
        }
        
        statflag = 0ull;
        if(b) {                                                                 // deterministic rule in first plane starting from gene-specified plane number
            statflag |= F_notgolrul;
            birthinitplane=golgij&0x3f;
            if(birthinitplane) b = (b>>birthinitplane) | (b<<(64-birthinitplane));    // rotate b to start at genetic birthplane
            FIRST1INDEX(b, kb);
            kb = (kb+birthinitplane)&0x3f;
            if((b3>>kb)&1ull) {                                                 // gene copied preferentially according to 2D birth rule with ancestor coming from
                for(k=0,nbmask=0ull;k<8;k++) {                                  // compute no of live neighbours using birthplane (coupled plane neighbours)
                    nbmask |= ((gol[deltaxy(ij,nb1x[k],nb1y[k])]>>kb)&1ull)<<k;
                }
                kch=selectdifft3(nbmask, &crot, &kodd);
            }
            else {                                                              // 3d-live-nb birth with gene copied from most difft member of neighbourhood
                plane = kb;                                                     // find plane with most live neighbours from the 3 planes kb + [0,1,-1]
                sum=(sums16[plane&0x3]>>(plane&0x3c))&0xf;
                if (sum<3) {                                                    // sum <3 in central plane, try next plane
                    plane = (kb+1) & 0x3f;
                    sum=(sums16[plane&0x3]>>(plane&0x3c))&0xf;
                }
                if (sum<3) {                                                    // sum <3 in two planes, try third
                    plane = (kb-1) & 0x3f;
                    sum=(sums16[plane&0x3]>>(plane&0x3c))&0xf;
                    if (sum<2) fprintf(stderr,"error in update_gol64 step %d ij %d, sum %d still <2 after 3 plane examination\n",totsteps,ij,sum);
                }
                for(k=0,nbmask=0ull;k<8;k++)                                    // compute live neighbour mask using chosen ancestor plane, excluding central site
                    nbmask |= ((gol[deltaxy(ij,nb1x[k],nb1y[k])]>>plane)&1ull)<<k;

                switch (sum-((gol[deltaxy(ij,nb1x[8],nb1y[8])]>>plane)&1ull)) { // do not take ancestor gene from central site on neighbour planes
                    case 1:  kch=selectdifft1(nbmask, &crot, &kodd); break;
                    case 2:  kch=selectdifft2(nbmask, &crot, &kodd); break;
                    case 3:  kch=selectdifft3(nbmask, &crot, &kodd); break;
                    case 4:  kch=selectdifft4(nbmask, &crot, &kodd); break;
                    case 5:  kch=selectdifft5(nbmask, &crot, &kodd); break;
                    case 6:  kch=selectdifft6(nbmask, &crot, &kodd); break;
                    default: kch=0;fprintf(stderr,"in update_gol64 error reached invalid case of sum %d\n",sum);
                }
            }
            ij1=deltaxy(ij,nb1x[kch],nb1y[kch]);
            newgene = golg[ij1];
            RAND128P(randnr);                                                   // inline exp so compiler recognizes auto-vec,
            rand2 = randprob(pmutmask, (unsigned int) randnr);
            nmut = (randnr >> 56) & 0x3f;                                       // choose mutation pos for length 64 gene : from bits 56:61 of randnr
            ancestor = newgene;
            newgene = newgene ^ (rand2<<nmut);                                  // introduce single mutation with probability pmut = probmut
            newgolg[ij]=newgene;
            statflag = statflag | F_birth;
            if (rand2) statflag = statflag | F_mutation;
            if(statflag&F_notgolrul) statflag |= F_nongolchg;
            if (golij) hashdeletegene(golgij,"step %d hash storage error 1 in update_gol64, gene %llx not stored\n"); // gene may be present because of live state on different plane
            hashaddgene(newgene,ancestor);
        } //end if birth
        else if (newgol[ij]) {                                                   // gene persistence if newgol non zero
            newgolg[ij]=golgij;
            statflag = statflag | F_survival;
        }
        else {                                                                   // gene death if was alive and death on all planes
            if(golij) {                                                          // site was alive, so it must have had a gene
                if(golgij==0ull) fprintf(stderr,"step %d ij %d golij non zero %llx but golgij zero\n",totsteps,ij,golij);
                hashdeletegene(golgij,"step %d hash storage error 2 in update_gol64, gene %llx not stored\n");
                newgolg[ij]=gene0;
                statflag = statflag | F_death;
            }
            else newgolg[ij]=golgij;                                           // stays dead
        } // end else
        newgolgstats[ij] = statflag;
        // fprintf(stderr,"step %d ij %d golij %llx b3 %llx b3d %llx golg[ij] %llx newgol[ij] %llx newgolg[ij] %llx\n",totsteps,ij,golij,b3,b3d,golg[ij],newgol[ij],newgolg[ij]);
    }
    
    if(randomsoup) random_soup(gol,golg,newgol,newgolg);
    if(vscrolling) v_scroll(newgol,newgolg);
    for (ij=0; ij<N2; ij++) {       // complete missing hash table records including activities, NB all sites have genes in gol16 but don't record activity for zerogene
        if(gol[ij]) hashgeneextinction(golg[ij],"hash storage error 4 in update_gol64, gene %llx not stored\n");
        if(newgolg[ij]!=gene0) hashgeneactivity(newgolg[ij],"hash storage error 5 in update_gol64, gene %llx not stored\n");
    }
    
    // fprintf(stderr,"step %d NbP %d planemask %llx\n",totsteps, NbP,  plcodingmask );
}
//------------------------------------------------------- update_gol2match ------------------------------------------------------------------------------
void update_gol2match(uint64_t gol[], uint64_t golg[],uint64_t newgol[], uint64_t newgolg[]) {     // selection model 14,15
// gol states on one plane interact with genes hard coupled to the gol states on a second plane
// this may also be viewed simply as genes with gol dynamics interacting with non-genetic entities with gol dynamics
// both gol planes are coded as two bottom bits in the gol[i,j] array
// coupling between the two planes occurs through complementary (1<->0) exact matching of the local Moore neighborhoods (8-sites)
// survival and/or birth rules are negated (opposite outcome) when coupling is established
    unsigned int ij,ij1,k,kch,nmut;
    uint64_t s0,s1,s0_2or3,s1_2or3,nbmask0,nbmask1,ancestor,newgene;
    uint64_t golmix,gol0,gol1,newgol1;
    uint64_t randnr, rand2, statflag;
    int nb[8],crot,kodd;
    int nb1x[8] = {-1,0,1,1,1,0,-1,-1};
    int nb1y[8] = {-1,-1,-1,0,1,1,1,0};
    
    canonical = repscheme & R_2_canonical_nb;
    
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
            if(s1&1ull) kch=selectdifft3(nbmask1, &crot, &kodd);                // 3 nbs
            else        kch=selectdifft2(nbmask1, &crot, &kodd);                // 2 nbs
            ij1=deltaxy(ij,nb1x[kch],nb1y[kch]);
            newgene = golg[ij1];
            RAND128P(randnr);                                                   // inline exp so compiler recognizes auto-vec,
            rand2 = randprob(pmutmask, (unsigned int) randnr);
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

    if(randomsoup) random_soup(gol,golg,newgol,newgolg);
    if(vscrolling) v_scroll(newgol,newgolg);
    for (ij=0; ij<N2; ij++) {       // complete missing hash table records including activities
        if(gol[ij]>>1) hashgeneextinction(golg[ij],"hash storage error 3 in update_gol2, gene %llx not stored\n");
        if(newgol[ij]>>1) hashgeneactivity(newgolg[ij],"hash storage error 4 in update_gol2, gene %llx not stored\n");
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
                case 0:
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
//.......................................................................................................................................................
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
    int nspecies,ngenealogydeep;
    int activitieshash(void);                                              // count activities of all currently active species
    int genealogies(void);                                                    // genealogies of all currently active species
    
    nhistG = nhist;
    nstatG = nstat;
    for (t=0; t<nsteps; t++) {
        newgol = planes[newPlane];
        newgolg = planesg[newPlane];
        newgolgstats = planesgs[newPlane];

        totsteps++;
        if(!(totsteps%10)) {
            nquadpatterns = hashtable_count(&quadtable);
            fprintf(stderr,"iteration step %d nquadpatterns %d\r",totsteps,nquadpatterns);
        }
        
        if (selection<8)        update(gol,golg,newgol,newgolg);              // calculate next iteration with selection
        else if (selection<10)  update_lut_sum(gol,golg,newgol,newgolg);      // calculate next iteration for lut sum version s= 1-8
        else if (selection<12)  update_lut_dist(gol,golg,newgol,newgolg);     // calculate next iteration for lut dist (corner/edge) version s=1-7
        else if (selection<14)  update_lut_canon_rot(gol,golg,newgol,newgolg);// calculate next iteration for lut canonical rotation version s=2-6
        else if (selection<16)  update_lut_2D_sym(gol,golg,newgol,newgolg);   // calculate next iteration for lut fully 2D symmetric version s=0-4
                                                                              // multiplane cases
        else if (selection<20) update_gol16(gol,golg,newgol,newgolg);         // calculate next iteration for 2-16x multiplane version : 1 or 2-4 other plane(s) coupled
        else if (selection<24) update_gol64(gol,golg,newgol,newgolg);         // calculate next iteration for 64x multiplane version
        else if (selection<25) update_gol2match(gol,golg,newgol,newgolg);     // calculate next iteration for 2xgol multiplane version

        if(nhist && (totsteps%nhist == 0)) countconfigs();                    // count configurations
        if(nstat && (totsteps%nstat == 0)) tracestats(gol,golg,golgstats,N2); // time trace point

        curPlane = (curPlane +1) % numPlane;                                  // update plane pointers to next cyclic position
        newPlane = (newPlane +1) % numPlane;
        gol = planes[curPlane];                                               // get planes of gol,golg data
        golg = planesg[curPlane];
        golgstats = planesgs[curPlane];

        nspecies=activitieshash();                                           // colors acttrace and returns current population arrays, need to run always for continuity
        if(nspecies<0) fprintf(stderr,"error returned from activitieshash\n");
        if(colorfunction>4 && colorfunction<8) {
            ngenealogydeep=genealogies();                                    // colors genealogytrace
            if(ngenealogydeep<0) fprintf(stderr,"error returned from genealogies\n");
        }
        
        totdisp++;                                                            // currently every step is counted for display in activities
    }
}
//------------------------------------------------------- initialize_planes ---------------------------------------------------------------------------
void initialize_planes(int offs[],  int No) {
    int i,j,idx;
    static int notfirst = 0;

    curPlane = 0;
    newPlane = 1;
    if (notfirst)   return;     // need to fix memory free at two levels unless this fix: no changes in planes structure during run allowed
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
    planesgs[0] = planegs0; planesgs[1] = planegs1;
#if maxPlane > 2
    planes[2]  = plane2;  planes[3]  = plane3;
    planesg[2] = planeg2; planesg[3] = planeg3;
    planesgs[2] = planegs2; planesgs[3] = planegs3;
#endif
#if maxPlane >4
    planes[4]  = plane4;  planes[5]  = plane5;  planes[6]  = plane6;  planes[7]  = plane7;
    planesg[4] = planeg4; planesg[5] = planeg5; planesg[6] = planeg6; planesg[7] = planeg7;
    planesgs[4] = planegs4; planesgs[5] = planegs5; planesgs[6] = planegs6; planesgs[7] = planegs7;
#endif
}
//------------------------------------------------------- readFile and writeFile ---------------------------------------------------------------
int readFile(char * code, const char *fileName) {
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
//.......................................................................................................................................................
int writeFile(char *fileName)  {   // initialize 32x32 genepat file with all empty sites
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
    unsigned int ncodingin;
    uint64_t g;
    static unsigned int rmask = (1 << 15) - 1;
    static int notfirst = 0;
    uint64_t startgenes[16];
    char *golgin;
    
    /* g=1ull;                                                              // test of FIRST1INDEX
    for(j=0;j<10;j++) {
        FIRST1INDEX(g,k);
        fprintf(stderr,"test of first1index cnt %d val %llx index %d\n",j,g,k);
        g*=42;
    } */

    /* g=1ull;                                                              // test of PATTERN4
    for(k=0;k<10;k++) {
        int found;
        for(j=0;j<16;j++) {
            PATTERN4(g,j,found);
            fprintf(stderr,"test of pattern4 pat %x val %llx found? %d\n",j,g,found);
        }
        g*=42;
    } */
    
    /* g=1ull;                                                              // test of PATTERN8
    for(k=0;k<10;k++) {
        int found;
        for(j=0;j<64;j++) {
            PATTERN8(g,j,found);
            fprintf(stderr,"test of pattern8 pat %x val %llx found? %d\n",j,g,found);
        }
        g*=42;
    } */
    
    srand(1234567); // Range: rand returns numbers in the range of [0, RAND_MAX ), and RAND_MAX is specified with a minimum value of 32,767. i.e. 15 bit
    state[0] = rand();state[1] = rand();
    cnt = 0;
    totsteps = 0;
    totdisp = 0;
    statcnts = 0;
    quadrants = -1;
    displayoneplane=64;
    rbackground = 0;
    nquadpatterns = 0;
    quadcollisions = 0;
    
    // writeFile("genepat.dat");                                            // can be used to initialize formatted template for gene input of 32x32 array

    rulemod = runparams[0];
    repscheme = runparams[1];
    selection = runparams[2];
    overwritemask = runparams[3];
    survivalmask = runparams[4];
    colorfunction = runparams[5];
    initfield = runparams[6];
    birthmask=runparams[7];
    
    randomsoup = 0;
    vscrolling = last_scrolled = 0;
    quadrants = -1;

    pmutmask = (unsigned int) simparams[0];                                 // low values of pmutmask <32 are interpreted as -nlog2pmut
    if(pmutmask<32) pmutmask = (0x1 << (32-pmutmask)) - 0x1ull;                  // NB if pmut==0, pmutmask==zero, no mutation.
    initial1density = simparams[1];
    initialrdensity = simparams[2];
    ncodingin = simparams[3];                                               // used in selection  4,5,6,8,
    ncoding = ncodingin & 0xff;
    ncoding2 = (ncodingin>>8) & 0xff;
    if (ncoding > 64) {
        fprintf(stderr,"value %d of ncoding is out of range\n",ncoding);
        ncoding = 64;
    }
    if (selection==8 || selection == 9) if (ncoding<1 || ncoding>4) ncoding = 4; // ncoding range restriction for selection = 8 ie 16*ncoding bits of gene used
    switch (selection) {
        case 16:
        case 17:
        case 18:
        case 19: NbP=(ncodingin >> 16)&0xf; if(NbP==0) NbP=16;NbG=NbP<<2;break;
        case 20:
        case 21:
        case 22:
        case 23: NbP=(ncodingin >> 16)&0x3f; if (NbP==0) NbP=64; NbG=NbP;break;
        case 24: NbP=2; NbG=64;break;
        default: NbP=1; NbG=64;
    }
    if (selection<8)
        codingmask = (1ull<<ncoding2)-1ull;                                 // bit mask corresponding to ncoding2 bits, only used in connection with add2ndmask1st
    else if (selection<10)
        codingmask = (1ull<<ncoding)-1ull;                                  // coding mask used to encode number of bits per LUT (1-4)
    startgenechoice = simparams[4];
    if (selection>=16 && selection<=19) displayplanes=0x1111111111111111ull;
    else if (selection>=20 && selection<24) displayplanes= ~0;
    
    fprintf(stderr,"___________________________________________________________________________________________\n");
    fprintf(stderr,"_________________________________ genelife simulation _____________________________________\n");
    fprintf(stderr,"runparams %d %d %d %d %d %d %d\n",runparams[0],runparams[1],runparams[2],
                                         runparams[3],runparams[4],runparams[5],runparams[6]);
    fprintf(stderr,"simparams %d %d %d %d %d\n",simparams[0],simparams[1],simparams[2],simparams[3],simparams[4]);
    fprintf(stderr,"pmutmask %x (NB 0 means no mutation) No of bit planes %d\n",pmutmask, NbP);

    if (selection == 16 || selection == 17) {                                // have to treat separately since bitplanes for gol only every 4th plane
        uint64_t pmask;
        pmask = (NbP==16) ? 0xffffffffffffffffull : (1ull<<(NbP<<2))-1;      // mask for bits included in np planes
        gene0=0xfedcba9876543210ull&pmask;                                   // default gene is to couple a plane only to itself
    }
    else                 gene0=0ull;                                         // normally default gene is 0 : unused when gol state not live

    if (selection >= 16 && selection <24) nstartgenes = 16;                 // need more start genes for these case but numbers arbitrary
    else nstartgenes = 8;

    switch (selection) {
        case 0:  for (k=0;k<4;k++) { startgenes[k] = 0xf0f0f0f0f0f0f0f0; startgenes[k+4] = 0x0f0f0f0f0f0f0f0f;} break;
        case 1:  for (k=0;k<8;k++)   startgenes[k] = ((0x1ull<<k*3)-1ull)<<20;break;
        case 2:
        case 3:  for (k=0;k<8;k++)   startgenes[k] =(((0x1ull<<20)-1ull)<<20)+((0x1ull<<k)-0x1ull);break;
        case 4:  for (k=0;k<8;k++)  {g = 0xff0ull; startgenes[k] = k<4 ? g+1 : ~((g<<16));} break;
        case 5:  for (k=0;k<8;k++)  {g = 0xf0ull + k; startgenes[k] = k<4 ? g : (~g)|(0xfull<<16);} break;
        case 6:
        case 7:  for (k=0;k<8;k++) startgenes[k]=(0x1ull<<(4+k*8))-1ull; break;
        
        case 8:  genegol[selection-8] = (codingmask<<((8+3-1)*ncoding))|(codingmask<<((2-1)*ncoding))|(codingmask<<((3-1)*ncoding));
                 for (k=0;k<8;k++)   startgenes[k] = genegol[selection-8];break;
        case 9:  genegol[selection-8] = 0xb32ull;                                 //GoL rule for survival in totalistic LUT case, variable length encoding
                 for (k=0;k<8;k++)   startgenes[k] = genegol[selection-8];  break;
        case 10: genegol[selection-8] = (0x7ull<<2)|(0xfull<<5)|(0xfull<<37);     //GoL rule for survival in corner/edge dist LUT case, fixed length encoding
                 for (k=0;k<8;k++)   startgenes[k] = genegol[selection-8];break;
        case 11: genegol[selection-8] = 0xbf3f27ull;                              //GoL rule for survival in corner/edge dist LUT case, variable length encoding
                 for (k=0;k<8;k++)   startgenes[k] = genegol[selection-8];break;
        case 12: genegol[selection-8] = 0xfull|(0x7full<<4)|(0x7full<<36);        //GoL rule for survival surv2 1st 4 surv3 next 7 birth3 7 bits from 32+4, canonical case
                 for (k=0;k<8;k++)   startgenes[k] = genegol[selection-8];break;
        case 13: genegol[selection-8] = 0x00cd00bbb73b3727;                       //GoL rule for S23 B3 with modular encoding for canonical rotation case
                 for (k=0;k<8;k++)   startgenes[k] = genegol[selection-8];break;
        case 14: genegol[selection-8] = (0x3full<<3)|(0x3ffull<<9)|(0x3ffull<<(32+9)); //GoL rule for S23 B3 for 2D_sym case with 1,2,6,10,13 configs for s=0,1,2,3,4
                 for (k=0;k<8;k++)   startgenes[k] = genegol[selection-8];break;
        case 15: genegol[selection-8] = 0x03f3ffb7b3373323;                       //GoL rule for S23 B3 with modular encoding for 2D_sym case
                 for (k=0;k<8;k++)   startgenes[k] =genegol[selection-8];break;
        case 16:
        case 17: for (k=0;k<16;k++)  startgenes[k] = gene0;                      // first set all startgenes as uncoupled
                 for (k=0;k<NbP;k++) startgenes[k] = (startgenes[k] & (~(0xfull<<(k<<2)))) | ((k+1)%NbP)<<(k<<2);break; // couple plane k+1 to k in startgene k
        case 18:
        case 19: for (k=0;k<16;k++)  startgenes[k] = gene0;                      // first set all startgenes as uncoupled
                 for (k=0;k<NbP;k++) startgenes[k] = 6<<(k<<2); break;           // coupled to neighbour plane before and after
        case 20:
        case 21: for (k=0;k<16;k++)  startgenes[k] = (k<<2)>NbP ? 0ull : (NbP==64) ? ~0ull : (1ull<<NbP)-1;break;
        case 24: for (k=0;k<8;k++)   startgenes[k] =((0x1ull<<k*3)-1ull)<<20;break;

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
    golgstats = planesgs[curPlane];

#ifdef HASH
    if(notfirst) {
        hashtable_term(&genetable);
        hashtable_term(&quadtable);
    }
    hashtable_init(&genetable,sizeof(genedata),N2<<2,0);     // initialize dictionary for genes
    hashtable_init(&quadtable,sizeof(quadnode),N2<<2,0);         // initialize dictionary for quadtree patterns
#endif
    notfirst = 1;
    if (initfield==1) {           // input from file genepat.dat with max size of 32*32 characters
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
    else {  // initfield gives linear size of random block for initialization (0 => full frame, as before)
        Nf = initfield;
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
                if(selection>=16 && selection<=19) for (k=0;k<NbP;k++) gol[ij] |= ((rand() & rmask) < initial1density)?(1ull<<(k<<2)):0ull;
                else for (k=0;k<NbP;k++) gol[ij] |= ((rand() & rmask) < initial1density)?(1ull<<k):0ull;
            }
        }
        for (ij=0; ij<N2; ij++) {
            g = 0ull;
            if (gol[ij]||((selection>=16)&&(selection<=19))) { // CHECK not 19 but 17? if live cell or multiplane, fill with random genome g or randomly chosen startgene depending on initialrdensity
                if (((unsigned) rand() & rmask) < initialrdensity) {for (k=0; k<NbG; k++) g = (g << 1) | (rand() & 0x1);g=gene0^g;}
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
    
#ifdef HASH
    for (ij=0; ij<N2; ij++) {
        if((gol[ij] && ((selection<16)||(selection>=20)))||((selection>=16)&&(selection<=19))||((gol[ij]>>1)&&(selection==24 || selection==25))) {
            hashaddgene(golg[ij],rootgene);
        }
    }
    // it is possible to enumerate keys and values
    hcnt=hashtable_count(&genetable);
    genotypes = hashtable_keys( &genetable );
    fprintf(stderr,"population size %d with %d different genes\n",cnt,hcnt);
    
    if(quadtree) qimage = quadimage(gol,N,0);
#endif
}
//------------------------------------------------------- set ...---------------------------------------------------------------------------
void set_colorfunction(int colorfunctionval) {
    if(colorfunction>9) fprintf(stderr,"error colorfunction value passed %d too large\n",colorfunctionval);
    else     colorfunction = colorfunctionval;
}
//.......................................................................................................................................................
int setget_act_ymax(int actymax) {                  // sets ymax for activities only if argument nonzero, reads old value
    int ymaxold;
    ymaxold = ymax;
    ymax = actymax;
    return(ymaxold);
}
//.......................................................................................................................................................
void set_selectedgene(uint64_t gene) {
    selectedgene=gene;
    fprintf(stderr,"selected gene set to %llx\n",selectedgene);
}
//.......................................................................................................................................................
void set_offsets(int dx,int dy,int dt) {
    offdx =dx;
    offdy = dy;
    if(dt>0) {
        dt=0;
        fprintf(stderr,"positive time offsets not allowed, looking into the future not possible\n");
    }
    if(dt<=-maxPlane) {
        dt=-maxPlane+1;
        fprintf(stderr,"not enough planes set for this time offset, recompile software with larger maxPlane value\n");
    }
    offdt = dt;
}
//.......................................................................................................................................................
void set_quadrant(int quadrant) {
    if (quadrant >= -1 && quadrant < 7) quadrants = quadrant;
    repscheme &= ~R_quadrant;                                           // remove all quadrant bits : only one set at a time in interactive version
    if(quadrant >= 0 && quadrant < 7) {
        repscheme |= R_14_quadrant_sele<<quadrant;                          // assumes quadrant selectors are 7 successive bits following R_14_...
    }
}
//.......................................................................................................................................................
void set_randomsoup(int randomsoupin) {
    randomsoup=randomsoupin;
}
//.......................................................................................................................................................
void set_rbackground(int rbackgroundin, int randomsoupin) {
    rbackground=rbackgroundin;
    randomsoup=randomsoupin;
    if(!rbackground) randomsoup=0;   // do not leave patch randomsoup variable active when turning off random background
}
//.......................................................................................................................................................
unsigned int set_repscheme_bits(int quadrant, int x, int y, int surviveover[]) {
    unsigned int quadrantval;
    
    quadrantval=(x<(Nmask>>1)? 0 : 1) + (y<(Nmask>>1)? 0 : 2);              // determine selected quadrant
    if(quadrants >= 0 && quadrants < 5) {                                   // assumes repscheme bits in pairs starting from bit 0 matching quadrants
        repscheme &=  ~(0x3llu<<(quadrants<<1));
        repscheme |=  ((uint64_t) quadrantval)<<(quadrant<<1);
    }
    else if (quadrant < 6) {
        survivalmask &= ~0x3u;
        survivalmask|= quadrantval;
    }
    else if (quadrant<7) {
        overwritemask &= ~0x3u;
        overwritemask|= quadrantval;
    }
    repscheme &= ~R_quadrant;                                                // remove all quadrant bits : only one set at a time in interactive version
    quadrants = -1;                                                          // reset internal quadrants choice so that full display is shown
    surviveover[0]=survivalmask;
    surviveover[1]=overwritemask;
    return(repscheme);
}
//.......................................................................................................................................................
void set_repscheme(unsigned int repscheme_in) {
    repscheme = repscheme_in;
}
//.......................................................................................................................................................
void set_rulemod(unsigned int rulemod_in) {
    rulemod = rulemod_in;
}
//.......................................................................................................................................................
void set_displayoneplane(unsigned int displayoneplane_in) {
    displayoneplane=displayoneplane_in;
}
//.......................................................................................................................................................
void set_displayplanes(unsigned int displayplanes_in) {
    int k;
    displayplanes=0ull;
    for(k=0;k<NbP;k++)
        displayplanes |= (uint64_t) (((displayplanes_in>>k)&0x1)<<(k<<2));
}
//.......................................................................................................................................................
void set_surviveover64(unsigned int surviveover[], int len ) {
    if (len==3) {
        survivalmask = surviveover[0];
        if(selection<8) overwritemask = surviveover[1];
        else {
            birthmask = surviveover[1];
            overwritemask = surviveover[2];
        }
    }
    else fprintf(stderr,"surviveover64 needs three parameters, %d provided\n",len);
}
//.......................................................................................................................................................
void set_vscrolling() {
    vscrolling=1-vscrolling;
    if(!vscrolling) last_scrolled = 0;
}
//------------------------------------------------------- get ... ---------------------------------------------------------------------------
//.......................................................................................................................................................
int get_log2N() {
    return(log2N);
}
//.......................................................................................................................................................
void get_curgol(uint64_t outgol[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
	    outgol[ij] = planes[curPlane][ij];
    }
}
//.......................................................................................................................................................
void get_curgolg(uint64_t outgolg[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
	    outgolg[ij] = planesg[curPlane][ij];
    }
}
//.......................................................................................................................................................
void get_acttrace(uint64_t outgolg[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolg[ij] = acttrace[ij];
    }
}
//.......................................................................................................................................................
void get_genealogytrace(uint64_t outgolg[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolg[ij] = genealogytrace[ij];
    }
}
//.......................................................................................................................................................
int get_nspecies() {
    int k,nspecies,nspeciesnow;
    nspecies = hashtable_count(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );

    for (k=0,nspeciesnow=0; k<nspecies; k++)
        if(geneitems[k].popcount) nspeciesnow++;
    return(nspeciesnow);
}
//.......................................................................................................................................................
void get_curgolgstats(uint64_t outgolgstats[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolgstats[ij] = golgstats[ij];                       // Note that golgstats is not dealt with in planes !
    }
}
//------------------------------------------------------- countspecies ---------------------------------------------------------------------------
int cmpfunc (const void * pa, const void * pb) {
   // return ( *(int*)pa - *(int*)pb );
   return ((*(const uint64_t *)pa > *(const uint64_t *)pb)  ? 1 : -1);
}

int cmpfunc1 ( const void *pa, const void *pb ) {
    const uint64_t *a = (const uint64_t *) pa;
    const uint64_t *b = (const uint64_t *) pb;
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
    fprintf(stderr,"%d\t%d\t\t%d\t\t%d\t\t%d\n",rulemod,repscheme,selection,overwritemask,survivalmask);
    fprintf(stderr,"pmutmask\tinit1\tinitr\tncoding\tstartchoice\n");
    fprintf(stderr,"%x\t\t%d\t%d\t%d\t%d\n",pmutmask,initial1density,initialrdensity,ncoding,startgenechoice);
}

void countspecies() {                                                       // counts current species without using hash tables
    countspecies1(gol, golg, N2);
}
//.......................................................................................................................................................
#ifdef HASH
int cmpfunc2 (const void * pa, const void * pb) {
   return ( genotypes[*(const int*)pa] > genotypes[*(const int*)pb] ? 1 : -1);
}

int cmpfunc3 (const void * pa, const void * pb) {
   return ( geneitems[*(const int*)pa].popcount < geneitems[*(const int*)pb].popcount ? 1 : -1);
}

void countspecieshash() {  /* counts numbers of all different species using qsort first */
    int k, *golgs, nspecies, nspeciesnow, nones;
    uint64_t last;
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

        if((genedataptr = (genedata *) hashtable_find(&genetable, last)) != NULL) {
            if(genedataptr->popcount) {
                fprintf(stderr,"count species %7d with gene %16llx has counts %7d and %3d ones\n",k,last,
                    genedataptr->popcount,nones);
                nspeciesnow++;
            }
        }
        else {
            fprintf(stderr,"countspecieshash popcount error, no entry in hash table\n");
            fprintf(stderr,"count species %d with gene %llx has counts ?? and %d ones\n",k,last,
                    nones);
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
int activitieshash() {  /* count activities of all currently active species */
    int i, j, ij, ij1, x, nspecies, nspeciesnow, popcnt0, popcnt1;
    int *gindices,*popln,*activities;
    double act;
    uint64_t *genes;
    const int maxact = 10000;
    // int ymax1;
    uint64_t gene;
    
    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );

    // if(gindices != NULL && col) free(gindices);
    for (i=0,nspeciesnow=0; i<nspecies; i++)
        nspeciesnow+=geneitems[i].popcount ? 1 : 0;
    
    gindices = (int *) malloc(nspeciesnow*sizeof(int));

    for (i=j=0; i<nspecies; i++) {
        if(geneitems[i].popcount) {
            gindices[j]=i;                                         // if col is 0 then the array gindices must be passed with sufficient length
            j++;
        }
        //gindices[j]=(geneitems[i].popcount) ? i : gindices[j--];   // if col is 0 then the array gindices must be passed with sufficient length
        //j++;
    }
    if (nspeciesnow > maxact) {                             //sort with in order of decreasing population
        qsort(gindices, nspeciesnow, sizeof(int), cmpfunc3);        // sort in decreasing count order
    }
    if (nspeciesnow > maxact) nspeciesnow = maxact;

    genes = (uint64_t *) malloc(nspeciesnow*sizeof(uint64_t));
    popln = (int *) malloc(nspeciesnow*sizeof(int));
    activities = (int *) malloc(nspeciesnow*sizeof(int));

    for (i=0; i<nspeciesnow; i++) {
        genes[i]=genotypes[gindices[i]];
        popln[i]=geneitems[gindices[i]].popcount;
        activities[i]=geneitems[gindices[i]].activity;
    }
    
    if (totdisp>=N) {                                               // 1 pixel to left scroll when full
        for(ij=0;ij<N2;ij++) {
            ij1 = ((ij+1)&Nmask)+((ij>>log2N)<<log2N);              // (i+1)%N+j*N;
            // if(ij1>=N2) fprintf(stderr,"error in scroll of acttrace\n");
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
    for(j=0;j<nspeciesnow;j++) {
        // activities[j] = N-1 - (activities[j] * (N-1)) / ymax;        // linear scale, needs truncation if ymax superceded
        act = (double) activities[j];
        // activities[j] = (N-1) - (int) ((N-1)*log2(act)/log2ymax);      // logarithmic scale, suffers from discrete steps at bottom
        activities[j] = (N-1) - (int) ((N-1)*act/(act+(double)ymax));
        gene = genes[j];
        ij = (x&Nmask)+activities[j]*N;
        // if(ij >= 0 && ij<N2)
        if(acttrace[ij]==rootgene)                                      // only one genotype to plot
            acttrace[ij] = gene;
        else {                                                          // plot species color with largest current population size if chpoice of multiple
            if((genedataptr = (genedata *) hashtable_find(&genetable, acttrace[ij])) != NULL)
               popcnt0 = genedataptr->popcount;
            else popcnt0 = 0;
            if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL)
               popcnt1 = genedataptr->popcount;
            else popcnt1 = 0;
            if(popcnt1 >= popcnt0) acttrace[ij]=gene;
        }
        //else
        //    fprintf(stderr,"error activity ij out of range activities[j] %d\n",activities[j]);
    }

    free(gindices);free(activities);free(genes);free(popln);
    return(nspeciesnow);
}
#else
int activitieshash() {return(0)}  /* dummy routine, no effect */
#endif

#ifdef HASH
int activitieshashx(int gindices[], uint64_t genes[], int popln[], int activities[]) {  /* count activities of all currently active species */
    int i, j, nspecies, nspeciesnow;
    const int maxact = 10000;
    
    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );

    // if(gindices != NULL && col) free(gindices);
    for (i=0,nspeciesnow=0; i<nspecies; i++)
        //if(geneitems[i].popcount) nspeciesnow++;                  // time critical loop, replace if by conditional
        nspeciesnow+= geneitems[i].popcount ? 1 : 0;

    if (nspeciesnow>10000) return(-1);                              // exit with error need to allocate more space in python
    for (i=j=0; i<nspecies; i++) {
        //if(geneitems[i].popcount) {
        //    gindices[j]=i;                                        // if col is 0 then the array gindices must be passed with sufficient length
        //    j++;
        //}
        gindices[j]=geneitems[i].popcount ? i : gindices[j--];      // time critical loop, replace if by conditional
        j++;
    }
    qsort(gindices, nspeciesnow, sizeof(int), cmpfunc3);            // sort in decreasing count order

    if (nspeciesnow > maxact) nspeciesnow = maxact;

    for (i=0; i<nspeciesnow; i++) {
        genes[i]=genotypes[gindices[i]];
        popln[i]=geneitems[gindices[i]].popcount;
        activities[i]=geneitems[gindices[i]].activity;
    }
    
    return(nspeciesnow);                                            // exit here unless doing display

}
#else
int activitieshashx(int gindices[], uint64_t genes[], int popln[], int activities[]) {return(-1)}  /* dummy routine, no effect */
#endif

int get_sorted_popln_act( int gindices[], uint64_t genes[], int popln[], int activities[]) {
    int nspecies;
    nspecies=activitieshashx(gindices, genes, popln, activities);        // sets acttrace and returns current population arrays
    return(nspecies);
}

//------------------------------------------------------- genealogies ---------------------------------------------------------------------------
#ifdef HASH
int cmpfunc4 (const void * pa, const void * pb) {
   return ( geneitems[*(const int *)pa].firstbirthframe > geneitems[*(const int *)pb].firstbirthframe ? 1 : -1);
}

int cmpfunc5 (const void * pa, const void * pb) {               // sort according to ancestry in genealogytrace

   int i1,i2,ij1,ij2,j;
   uint64_t gene1,gene2;
   i1=*(const int *)pa; i2=*(const int *)pb;
   
   for (j=0;j<genealogydepth;j++) {
        ij1 = i1+j*N; ij2 = i2+j*N;
        gene1=working[ij1]; gene2=working[ij2];
        if(gene1!=gene2) return((((gene1 > gene2) && (gene1!=rootgene)) || (gene2==rootgene)) ? 1 : -1);
    }
    return(0);
}

int cmpfunc6 (const void * pa, const void * pb) {               // sort according to ancestry in genealogytrace using activity ordering
   int i1,i2,ij1,ij2,j;
   uint64_t gene1,gene2;
   i1=*(const int *)pa; i2=*(const int *)pb;
   int act1,act2;
   
   for (j=0;j<genealogydepth;j++) {
        ij1 = i1+j*N; ij2 = i2+j*N;
        gene1=working[ij1]; gene2=working[ij2];
        if(gene1!=gene2)  {
            if((gene1!=rootgene) && (gene2 != rootgene)) {
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene1)) != NULL) act1 = genedataptr->activity; else act1 = 0;
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene2)) != NULL) act2 = genedataptr->activity; else act2 = 0;
                return(act1 > act2 ? 1 : (act1==act2 ? (gene1 > gene2 ? 1 : -1) : -1));
            }
            else return((gene2==rootgene) ? 1 : -1);
        }
    }
    return(0);
}

int cmpfunc7 (const void * pa, const void * pb) {               // sort according to ancestry in genealogytrace using popln size ordering
   int i1,i2,ij1,ij2,j;
   uint64_t gene1,gene2;
   i1=*(const int *)pa; i2=*(const int *)pb;
   int pop1,pop2;
   
   for (j=0;j<genealogydepth;j++) {
        ij1 = i1+j*N; ij2 = i2+j*N;
        gene1=working[ij1]; gene2=working[ij2];
        if(gene1!=gene2)  {
            if((gene1!=rootgene) && (gene2 != rootgene)) {
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene1)) != NULL) pop1 = genedataptr->popcount; else pop1 = 0;
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene2)) != NULL) pop2 = genedataptr->popcount; else pop2 = 0;
                return(pop1 > pop2 ? 1 : (pop1==pop2 ? (gene1 > gene2 ? 1 : -1) : -1));
            }
            else return((gene2==rootgene) ? 1 : -1);
        }
    }
    return(0);
}

int genealogies() {  /* genealogies of all currently active species */
    int j, jmax, i, ij, nspecies, nspeciesnow, birthstep;
    int j1, j2, j3, activity, gorder[N];
    uint64_t gene, ancgene, nextgene;
    int *gindices,*popln,*activities,*birthsteps;
    uint64_t *genes;

    
    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );

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
    birthsteps = (int *) malloc(nspeciesnow*sizeof(int));

    for (i=0; i<nspeciesnow; i++) {
        genes[i]=genotypes[gindices[i]];
        popln[i]=geneitems[gindices[i]].popcount;
        activities[i]=geneitems[gindices[i]].activity;
        birthsteps[i]=geneitems[gindices[i]].firstbirthframe;
    }
    
    for(ij=0;ij<N2;ij++) working[ij]=rootgene;              // set field to rootgene black
    ancgene=rootgene;                                             // never really used, but included to avoid unitialized warning
    birthstep=0;
    activitymax=0;
    for (i=jmax=0; i<nspeciesnow; i++) {
        //j1=0;
        for (j=0;j<N;j++) {  // go back at most N links in genealogy
            if(j) {
                gene=ancgene;
                if(gene==rootgene) {
                    ancgene = rootgene;
                    // birthstep = 0;  // not needed
                    // activity = 0;   // not needed as does not affect max value
                }
                else {
                    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
                        ancgene=genedataptr->firstancestor;
                        // birthstep = genedataptr->firstbirthframe; // not needed
                        activity = genedataptr->activity;
                        if(activity>activitymax) activitymax=activity;
                    }
                    else fprintf(stderr,"ancestor not found in genealogies\n");
                }
            }
            else  {
                gene=genes[i];
                ancgene=geneitems[gindices[i]].firstancestor;
                // birthstep=geneitems[gindices[i]].firstbirthframe; // not needed
                activity=geneitems[gindices[i]].activity;
                if(activity>activitymax) activitymax=activity;
            }
            if (gene == rootgene) break;                            // reached root, exit j loop
            ij = i+j*N;
            working[ij]=gene;
        }
        if (j>jmax) jmax=j;
    }
    genealogydepth = jmax;
    
                                                                    //reverse ancestries to allow comparison at same number of speciations
    for (i=0; i<nspeciesnow; i++) {
        for(j=0;j<N;j++) {
            if (working[i+j*N]==rootgene) break;
        }
        for(j1=0;j1<(j>>1);j1++) {
            gene=working[i+(j-j1-1)*N];
            working[i+(j-j1-1)*N]=working[i+j1*N];
            working[i+j1*N]=gene;
        }
    }
    for (i=0; i<N; i++) gorder[i]=i;
    qsort(gorder, nspeciesnow, sizeof(int), cmpfunc5);              // sort according to ancestral lines - use cmpfunc5 to sorting genes laterally via gene value
    //qsort(gorder, nspeciesnow, sizeof(int), cmpfunc6);            // sort according to ancestral lines - use cmpfunc6 to sorting genes laterally via activity
    //qsort(gorder, nspeciesnow, sizeof(int), cmpfunc7);            // sort according to ancestral lines - use cmpfunc7 to sort genes laterally via population size
    
    for (i=0;i<N;i++) if((gorder[i]<0)||(gorder[i]>=N)) fprintf(stderr,"step %d error in gorder out of bounds at i = %d with value %d\n",totsteps,i,gorder[i]);

    for(ij=0;ij<N2;ij++) genealogytrace[ij]=rootgene;               // initialize genealogytrace to root gene before drawing part of it

    if((colorfunction==6)||(colorfunction==7)) {                                          // time trace of genealogies
      for(i=0;i<nspeciesnow;i++) {
        for(j=0,j1=0;j<jmax;j++) {
            if(gorder[i]>=nspeciesnow) fprintf(stderr,"error in genealogies gorder at i=%d, order value %d out of range\n",i,gorder[i]);
            ij = gorder[i]+j*N;
            gene = working[ij];
            ij+=N;
            if(ij<N2) {
                nextgene = working[ij];
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
      // for(i=nspeciesnow;i<N;i++) for(j=0;j<N;j++) genealogytrace[gorder[i]+j*N]=rootgene;
    }
    else {                                                          // species changes only trace (colorfunction == 5)
      for(i=0;i<nspeciesnow;i++) {
        for(j=0;j<jmax;j++) {
            ij=i+j*N;
            genealogytrace[ij]=working[gorder[i]+j*N];
        }
      }
    }
    free(gindices);free(activities);free(genes);free(popln);free(birthsteps);

    return(jmax);
}
#else
int genealogies(int gindices[], uint64_t genes[], int popln[], int activities[], int birthsteps[]) {return(0)}  /* dummy routine, no effect */
#endif

//------------------------------------------------------- misc ---------------------------------------------------------------------------
void delay(int milliseconds) {
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
