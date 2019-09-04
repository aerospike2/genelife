// subgenelife.c
// Written by John S. McCaskill and Norman H. Packard
//
// Project fastgenegol
//
// First created by John McCaskill on 14.07.17. Last modified Aug 2019.
// Copyright © 2017,2018,2019 European Center for Living Technology. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>

//-----------------------------------------------------------size of array -------------------------------------------------------------------------
const int log2N = 9;                // toroidal array of side length N = 2 to the power of log2N (minimum log2N is 6 i.e. 64x64)
const int N = 0x1 << log2N;         // only side lengths powers of 2 allowed to enable efficient implementation of periodic boundaries
const int N2 = N*N;                 // number of sites in square-toroidal array
const int Nmask = N - 1;            // bit mask for side length, used instead of modulo operation
const int N2mask = N2 - 1;          // bit mask for array, used instead of modulo operation
const uint64_t rootclone = N2;      // single bit for mask enquiries for clones with root heritage
//--------------------------------------------------------- main parameters of model ----------------------------------------------------------------
unsigned int rulemod = 1;           // determine whether to modify GoL rules
int selection = 1;                  // 0-7 fitness of 2 live neighbours, 8-15 lut symmetry and encoding scheme
                                    // 0.  integer value                        1. number of ones         2. scissors-stone-well-paper: wins over left 1-1-2-2
                                    // 3.  scissors-stone-well-paper 1-1-1-1    4. 5. predator prey       6. two target coding     7. selection always fails (no result)
                                    // 8.  genes encode lut based on sum s for survival/birth:  s 1-8 (8 configs): ncoding determines nr of fixed position bits per LUT entry
                                    // 9.  like 8 but position-independent encoding of lut entries: ncoding 1 bit only
                                    // 10. genes encode lut based on sum and edge-centred count (s,se) for survival/birth: s 1-7 se 0-4 (2+3+4+5+4+3+2=23 configs): ncoding 1-bit only
                                    // 11. like 10 but position-independent encoding of lut entries
                                    // 12. genes encode lut based on sum s and canonical rotations survival/birth sum 2-6 (4,7,10,7,4 canon.rotns) : ncoding 1 bit only
                                    // 13. like 12 but but position-independent encoding of lut entries
                                    // 14. genes encode lut based on 2D rotation/reflection symmetries for sum 0-4 only (1,2,6,10,13 cann.rotns) 2x32 bits: ncoding 1 bit only
                                    // 15. like 14 but position-independent encoding of lut entries
unsigned int repscheme = 1;         // replication scheme with separate meaning for selection 0-7 (update_23) and 8-15 other update fns
                                    // for selection 0-7: see #define R_... 1st section below
                                    // --  lowest 10 bits define 5 pairs of bit options for:
                                    // --  0,1 select birth 2,3 neighbour choice 4,5 enforce birth 6,7 2nd,1st neighbour genes 8,9 no successive nonGoL
                                    // --  bits 10-13 specify 4-bit mask for 2-live neighbour birth onl for those of the 4 canonical configurations set
                                    // --  bits 14-20 activate quadrant exploration of specific subset of the first 5 pairs of repscheme and survival and overwrite masks
                                    // --  bit 21 parentdies
                                    // for selection 8-15: see #define R... 2nd section below
unsigned int survivalmask = 0x3;    // for selection 0-7 survive mask for two (bit 1) and three (bit 0) live neighbours
                                    // for selection=8,9 it is 16-bit and for 10-11 it is 19-bit and 12-13 32-bit survival mask to restrict luts
unsigned int birthmask = 0x0;       // for selection 0-7 unused, birth control is done via bits 0,1,4,5 of repscheme
                                    // for selection=8,9 it is 16-bit and for 10-11 it is 19-bit and 12-13 32-bit birthmask to restrict luts
unsigned int overwritemask = 0x3;   // for selection 0-7 bit mask for 4 cases of overwrite: bit 0. s==3  bit 1. special birth s==2
                                    // for selection=8,9 it is 16-bit and for 10-11 it is 19-bit and 12-13 32-bit birthmask to restrict luts
unsigned int ancselectmask = 0xff;  // whether to use selection between genes to determine ancestor: for each birth rule (0 use positional choice via selectdifft)
int ncoding = 1;                    // byte 0 of python ncoding : number of coding bits per gene function
int ncoding2 = 0;                   // byte 1 of python ncoding: number of coding bits per gene function for masks in connection with repscheme add2ndmask1st R_6,7
unsigned int pmutmask;              // binary mask so that prob of choosing zero is pmut = pmutmask/2^32. If value<32 interpret as integer -log2(prob)
//...........................................................diagnostic control..........................................................................
const unsigned int diag_all                 = 0xffff;   // all diagnostics active
//const unsigned int diag_all               = 0xffdb;   // all diagnostics active except clones
unsigned int diagnostics                    = diag_all; // bit mask for diagnostics as defined by following constants
const unsigned int diag_hash_genes          = 0x1;      // enable hash storage of all genes encountered in simulation
const unsigned int diag_hash_patterns       = 0x2;      // enable hash storage of all patterns encountered in simulation
const unsigned int diag_hash_clones         = 0x4;      // enable hash storage of all clones encountered in simulation
const unsigned int diag_activities          = 0x8;      // enable activity statistics recording for genes and patterns
const unsigned int diag_simp_genealogies    = 0x10;     // enable simple genealogies with first or most recent ancestor
const unsigned int diag_clone_genealogies   = 0x20;     // enable genealogies by hashed clones
const unsigned int diag_component_labels    = 0x40;     // enable spatial connected component labelling, mapping and coloring
const unsigned int diag_offset_statistics   = 0x80;     // enable collection of offset histogram statistics: histo,numHisto,offsets,Noff
const unsigned int diag_scrolling_trace     = 0x100;    // enable scrolling time tracing of activities for genes and patterns,populations,genealogies
const unsigned int diag_longtime_trace      = 0x200;    // enable longer time tracing of activities for genes and patterns,populations (poss.genealogies)
const unsigned int diag_general_statistics  = 0x400;    // enable collection of general statistics: livesites,genestats,stepstats,configstats
const unsigned int diag_info_transfer_hist  = 0x800;    // enable collection of glider information transfer histogram in 8 directions  N E S W NE SE SW NW
//-----------------------------------------------------------initialization and color parameters---------------------------------------------------------
int initial1density = (1<<15)>>1;   // initial density of ones in gol as integer value, divide by 2^15 for true density
int initialrdensity = (1<<15)>>1;   // initial density of random genes in live sites, divide by 2^15 for true density
int startgenechoice = 8;            // selection for defined starting genes 0-8 (8 is random 0-7) otherwise choose only particular startgene
int initfield = 0;                  // 0 input from random field of random or start genes, 1 input from file genepat.dat of 32x32 indexes to start genes
                                    // value i>1: initialized to random field on central ixi block only, outside this zero.
                                    // value<0: initialize with gol and golg values
int colorfunction = 0;              // color function choice of 0: hash or 1: functional (color classes depends on selection parameter)
                                    // 2: as in 1 but color sites where last step was non GoL rule yellow, 3: as in 2 but yellow if state produced by non GoL
                                    // 4: activities 5: populations 6: genealogies without time 7: genealogies with time brightness by activity 8: gliders
                                    // 9: connected components 10 : connected component activities 11: genealogy based individual colors 12: genetic glider det
int colorfunction2 = -1;            // colorfunction for second window: as above, but -1 means same value as colorfunction : only one fn call made
int colorupdate1 = 1;               // flag to enable routine print statements to terminal during run : linked to colorfunction display in python
int connectedprints = 0;            // flag to enable printouts of statistics of connected components (every 10 steps)
int ancestortype = 0;               // display and return genealogies via first ancestor (0), clonal ancestry (1) or first & clonal in 2 win (2)
int parentdies = 0;                 // model variant enhancing interpretation of non-proliferative birth as movement (1) or default (0): set in repscheme
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
#define R_14_parentdies_23 0x4000   /* 1. parent is forced to die on birth 0: not: only used for selection 0-7 ie update_23 */
#define R_15_dummy         0x8000   /* not yet used but connected up to graphical interface */
#define R_quadrant         0x7f0000 /* 1: quarter the spatial domain with one or more of 7 pairs of repscheme bits ie 4 different values */
#define R_16_quadrant_sele 0x10000  /* q0 1: quarter the spatial domain with selection enable values for 2,3 live nbs: only in update ie for selection<8 */
#define R_17_quadrant_posn 0x20000  /* q1 1: quarter the spatial domain with selection enable values for 2,3 live nbs: only in update ie for selection<8 */
#define R_18_quadrant_enfb 0x40000  /* q2 1: quarter the spatial domain with enforce birth values for 2,3 live nbs: only in update ie for selection<8 */
#define R_19_quadrant_2nb1 0x80000  /* q3 1: quarter the spatial domain with 1st nb masks and/or 2nd nb addition: only in update ie for selection<8 */
#define R_20_quadrant_ngol 0x100000 /* q4 1: quarter the spatial domain with last non gol rule and/or non gol created state: only in update ie for selection<8 */
#define R_21_quadrant_surv 0x200000 /* q5 1: quarter the spatial domain with survival values for 2,3 live nbs: only in update ie for selection<8 */
#define R_22_quadrant_over 0x400000 /* q6 1: quarter the spatial domain with overwrite values for 2,3 live nbs: only in update ie for selection<8 */
//.................................................. LUT repscheme (selection 8-15) repscheme bits .......................................................
#define R_0_survivalgene  0x1       /* 1: survival gene chosen from central existing gene 0: survival gene taken from neighbours as in birth */
#define R_1_nb_OR_AND     0x2       /* 1: OR of neighbours determines genetic LUT in selection 8,10,12,14 0: AND of neighbours */
#define R_2_canonical_nb  0x4       /* 1: choose live neighbour at zero bit in canonical rotation 0: choose most difft position */
#define R_3_parentdies    0x8       /* 1: parent is forced to die on birth 0 not. Only used for selection 8-15 */
#define R_46_repselect    0x70      /* 0-7 choice of selection mechanism for LUT genes : 0: min 1: max 2: min 1s 3: max 1s 4: neutral 5: neutral difft 6,7: c-S-2B */
#define R_47_repselect    0xf0      /* +8-15 choice of golr selection mechanism : 8: 9: 10: 11: 12: 13: 14: 15:  */
#define R_11_golr_select  0x80      /* 1: choose 8-15 above 0: choose 0-7 above */
#define R_810_disambig    0x700     /* 0-7 choice of different disambiguation mechanisms for symmetric canonical rotations */
#define R_7_random_resln  0x800     /* 1: random choice amongst selected live neighbours 0: deterministic choice based on gene content and position */
//----------------------------------------status flag bits for recording site status in golgstats array---------------------------------------------------
const int F_s_live =      0x7ull;      // s value mod 8 (number of live neighbors) for selection 8-15  and separate bits below for selection 0 to 7
const int F_1_live =      0x1ull;      // bit is bit0 of s for selection 8-15 or 1 if exactly 1 live neighbours for selection 0-7 : currently not set
const int F_2_live =      0x2ull;      // bit is bit1 of s for selection 8-15 or 1 if exactly 2 live neighbours for selection 0-7
const int F_3_live =      0x4ull;      // bit is bit2 of s for selection 8-15 or 1 if exactly 3 live neighbours
const int F_birth =       0x8ull;      // bit is 1 if birth (includes overwriting of genes for some parameter values)
const int F_mutation =    0x10ull;     // bit is 1 if a mutation event occured
const int F_disambig =    0x80ull;     // bit is 1 if the disambiguate routine is used: nbest>1, nsame > 1
const int F_survival =    0x20ull;     // bit is 1 if last step was a 1->1 gol survival
const int F_death =       0x40ull;     // bit is 1 if last step involved the death of a gene ie 1->0 gol transition
const int F_golstate =    0x100ull;    // bit is 1 if gol state is 1
const int F_golchange =   0x200ull;    // bit is 1 if state changed at last step
const int F_nongolchg =   0x400ull;    // bit is 1 if state when produced (ie changed to) was made by a non GoL rule
const int F_notgolrul =   0x800ull;    // bit is 1 if last step not a GoL rule
const int F_survmut =     0x1000ull;   // bit is 1 if mutation or survival from non-replicated mutant: mutation(t) or mutation(t-1)&survival(t)
const int F_parent =      0x2000ull;   // bit is 1 if individual that was at this site was parent/ancestor/genetic donor of a new individual born in last step
const int F_parentaldeath=0x4000ull;   // bit is 1 if individual died at this step only because it was a parent under parentdies=1 option
const int F_dummy =       0x8000ull;   // free: not yet used
const int F_livenbs =     0xff0000ull; // mask for storing configuration of live neighbours : clockwise from top-left neighbour (NW) (sel 0-7 for s>1, sel 8-15 for s>0)
const int F_2select =     0x1000000ull;// bit is 1 if the 2 live neighbour selection routine was employed : selection 0-7 only
const int F_3g_same =     0x2000000ull;// bit is 1 if exactly 3 live nbs and all 3 have same gene : only for selection 0-7
//................................................................ miscellaneous ......................................................................
#define ASCII_ESC         27        /* escape for printing terminal commands, such as cursor repositioning : only used in non-graphic version */
#define IJDEBUG       105690        /* ij nr for debug printouts */
//----------------------------------------------------------hash table implementation of python style dictionary---------------------------------------
#define HASHTABLE_IMPLEMENTATION    // uses Mattias Gustavsson's hashtable (github) for unsigned 64 bit key dictionary
#define HASHTABLE_U64 uint64_t      // define the hashtable 64bit unsigned int type as uint64_t
#define HASHTABLE_U32 uint32_t      // define the hashtable 32bit unsigned int type as uint32_t
#define HASHTABLE_SIZE_T uint64_t   // use 64 bit unsigned key type consistent with this file
#include "hashtable.h"              // Gustavsson's file was modified because the 64bit to 32bit key compression code produces 0 for some values
hashtable_t genetable;              // genetable contains the hash table for genes : the hash table keys are the 64bit genes themselves
typedef struct genedata {           // value of keys stored for each gene encountered in simulation
    unsigned int popcount;          // initialized to 1
    short unsigned int firsttime;   // first time the gene was created from ancestor: initialized to 0
    short unsigned int recenttime;  // most recent time the gene was created by mutation from ancestor or randomly generated : initialized to 0
    short unsigned int lasttime;    // last time this gene was seen : enables genes to be mapped to their last epoch (including present)
    short int lastextinctiontime;   // this is initialized to -1, meaning no extinctions yet
    unsigned int activity;          // initialized to 0
    unsigned int nextinctions;      // initialized to 0
    unsigned int dummy;             // pad to 64 bit word boundary
    uint64_t gene;                  // stored gene : note that two difft 64-bit genes may be stored at same location, so need check
    uint64_t firstancestor;         // this is initialized to a special gene seq (rootgene) not likely ever to occur for starting genes
} genedata;
const uint64_t rootgene = 0xfedcba9876543210; // initial special gene as root for genealogies
const uint64_t generepeat = 0x0123456789abcdef; // special gene sequence to denote repeated genes found in genealogy
genedata ginitdata = {1,0,0,0,-1,0,0,0,0ull,rootgene};  // initialization data structure for gene data
genedata *genedataptr;              // pointer to a genedata instance
HASHTABLE_SIZE_T const* genotypes;  // pointer to stored hash table keys (which are the genotypes)
genedata* geneitems;                // list of genedata structured items stored in hash table
int genefnindices[1<<24];           // table of activities for functional gene indices calculated by genefnindex
//.......................................................................................................................................................
hashtable_t quadtable;              // hash table for quad tree of spatial live state patterns
typedef struct quadnode {           // stored quadtree binary pattern nodes for population over time (currently only for analysis not computation)
    uint64_t hashkey;               // hash table look up key for node : enables tree exploration more directly than construction from nw,ne,sw,se including collision avoidance
//    uint64_t firstancestor;         // pattern immediately preceding the generation of this pattern
    uint64_t nw, ne, sw, se;        // constant keys to hashed quadnodes or 64-bit patterns : we terminate one level higher than Gosper & golly
    unsigned short int isnode;      // 1 if this is a node not a pattern
    unsigned short int size;        // side length of square image corresponding to quadtree pattern
    unsigned int activity;          // number of references to finding this node
    unsigned int pop1s;             // 32 bit number of 1s in quadnode
    unsigned int firsttime;         // first time node was identified
    unsigned int lasttime;          // last time node was identified : used to determine if part of current timestep
    unsigned int topactivity;       // activity of pattern as top of connected component (last field : giving size of record an even number of 64bit words)
} quadnode;
quadnode quadinit = {0ull,/*rootgene,*/0ull,0ull,0ull,0ull,0,0,1,0,0,0,0};
uint64_t qimage;                    // key for quadimage datastructure of entire image
int quadcollisions = 0;
HASHTABLE_SIZE_T const* quadkeys;   // pointer to stored hash table keys (which are the quadkeys)
quadnode* quaditems;                // list of quadnode structured items stored in hash table
typedef struct smallpatt {          // stored binary patterns for 4*4 subarrays or smaller (16bit)
    // unsigned short int size;        // side length of square image corresponding to pattern
    unsigned int topactivity;       // number of references to finding this pattern as top level of connected component
    unsigned int activity;          // number of references to finding this pattern
    unsigned int firsttime;         // first time pattern was identified
    unsigned int lasttime;          // last time pattern was identified
} smallpatt;
smallpatt smallpatts[65536];
//.......................................................................................................................................................
hashtable_t clonetable;             // hash table for clone ancestry
typedef struct clonedata {
    uint64_t birthid;               // birth time and place of clone (upper 32 bits: time, lower 32 bits: place ij) plus 1 bit (N2) 0: ancesotr 1: no ancestor
    uint64_t parentid;              // ancestor's id : i.e. birth time and place (upper 32 bits: time, lower 32 bits: place ij)
    uint64_t gene;                  // gene of this individual
    unsigned int popln;             // number of live individuals in this clone (32 bit)
    unsigned int activity;          // activity of clone (32 bit)
    // unsigned int subclones;      // number of mutant clones stemming from clone (32 bit): not used currently
} clonedata;
clonedata cinitdata = {0ull,0ull,0ull,1,1}; // initial values of clonedata: popln set to 1
clonedata *clonedataptr;            // pointer to clonedata instance
HASHTABLE_SIZE_T const* clones;     // pointer to stored hash table keys (which are the clones, live or with live descendants)
clonedata* cloneitems;              // list of clonedata structured items stored in hash table
//-------------------------------------------------------------------------------------------------------------------------------------------------------
int totsteps=0;                     // total number of simulation steps
int totdisp=0;                      // total number of displayed steps
int statcnts=0;                     // total number of statistic timepts counted
uint64_t codingmask;                // ncoding derived mask for ncoding bits
int nhistG = 0;                     // interval for collecting config histogram data : 0 no collection, nstatG collection with time
int nstatG = 0;                     // interval for collecting other statistical trace data : 0 no collection
int genealogydepth = 0;             // depth of genealogies in current population
int genealogycoldepth = 0;          // genes coloured by colour of ancestor at this depth in colorfunction=11
int ngenealogydeep;                 // depth of genealogy
int clonealogydepth = 0;            // depth of clonealogies in current population
int nclonealogydeep;                // depth of clonealogy
int ambigsum;                       // temporary statistic of ambiguous resolution cases
int nbshist[8];                     // global histogram to count neighbourhoods to check for assymetries
//.........................................................resource management NYI.......................................................................
int rmax = 1;                       // max number of resources per cell : 0 also turns off resource processing
int rthresh = 3;                    // minimum resource number per neighborhood to allow birth process
int rbirth = 2;                     // minimum resource number per neighborhood to allow birth replenishment
//---------------------------------------------------------main arrays in simulation---------------------------------------------------------------------
uint64_t *gol, *golg, *golb, *golr; // pointers to one plane of gol, golg, golb & golr arrays: live/dead, gene, cloneid (birth t,x,y), resource
uint64_t *golgstats;                // pointer to 64 bit masks for different events during processing at one of the plane cycle locations
uint64_t stashgol[N2];              // for stashing state and recovering with initfield = -1
uint64_t stashgolg[N2];             // for stashing state and recovering with initfield = -1
uint64_t golmix[N2];                // array for packing configs
uint64_t gene0;                     // uncoupled planes background gene, non zero for selection==16,17
uint64_t selectedgene;              // gene currently selected interactively in graphics window
uint64_t genegol[16];               // genes encoding for GoL in various LUT selection models indexed as selection-8
//---------------------------------------------------------additional variables set in simulation--------------------------------------------------------
unsigned int canonical;             // current value of choice of canonical position repscheme bit 2 : needed globally in ...difft2-6 routines
int quadrants=-1;                   // integer choice of bit pair from repscheme/survivalmask/overwritemask for quadrant division of array (-1 none)
int randominflux=0;                 // 1,2 steady or intermittent random input of genes and gol states into the square central region defined by initfield
                                    // 3 steady deletion of genes at random rate rbackground
int rbackground=0;                  // integer background rate of random gene input per frame per site : rate is rbackground/32768 (nonzero overides influx)
                                    // the gene input depends on randominflux value: 2 GoL 1 random
int vscrolling=0;                   // whether to do vertical scrolling to track upwards growth (losing all states that fall off downward cliff)
int vscrolly=0;                     // cumulative extent of scrolling (mod N) : to allow clone birthids to be parsed correctly
int last_scrolled = 0;              // whether vscrolling applied on last time step (needed for correct glider detection)
int ymax = 2000;                    // gene activity scale max for plotting : will be adjusted dynamically or by keys
int ymaxq = 2000;                   // quad pattern activity scale max for plotting : will be adjusted dynamically or by keys
// double log2ymax = 25.0;          // activity scale max 2^25 = 33.5 * 10^6 : suffers from discrete steps at bottom, not used
int activitymax;                    // max of activity in genealogical record of current population
int activityfnlut = 0;              // whether to lump activities of genes which have the same sequence at functionally masked in active positions of LUT
int noveltyfilter = 0;              // novelty filter for colorfunction 9 : if on (key "n"), darkens non-novel components (activity>1) in display
int activity_size_colormode = 0;    // color by size for colorfunction 10 : if on (key "p")  1 log2 enclosing square size 2 use #pixels 3 use sqrt(#pixels)
int xdisplay,ydisplay = -1;         // display x and y coordinates selected by mouse in python
int info_transfer_h = 0;            // whether to display histogram on glider information transfer counts (non zero)
int it_nbhood = 7;                  // size of neighborhood for collecting glider characterization histogram : default 7x7 nbhood
uint64_t gliderinfo[408];           // histogram of counts for glider detection by match quality in eight directions N E S W NE SE SW NW
//------------------------------------------------ arrays for time tracing, activity and genealogies ----------------------------------------------------
const int startarraysize = 1024;    // starting array size (used when initializing second run)
int arraysize = startarraysize;     // size of trace array (grows dynamically)
int *livesites = NULL;              // dynamic array pointer for statistics of number of live sites over time
int *genestats = NULL;              // dynamic array pointer for statistics of number of 4 genotype classes over time
int *stepstats = NULL;              // dynamic array pointer for statistics of site update types over time
int *configstats = NULL;            // dynamic array pointer for statistics of gol site configurations (x,y,t) offsets
uint64_t poptrace[N2];              // scrolled trace of last N time points of population of N most frequent genes
uint64_t acttrace[N2];              // scrolled trace of last N time points of activity of N most frequent genes
uint64_t acttraceq[N2];             // scrolled trace of last N time points of activity of N most frequent quad patterns
unsigned char acttraceqt[N2];       // type of entry in acttraceq : 1 quad, 0 smallpatt (<65536) i.e. corresponding to isnode
uint64_t genealogytrace[N2];        // image trace of genealogies for N most frequently populated genes
uint64_t clonealogytrace[N2];       // image trace of clonealogies for N most frequently populated clones
const int nNhist = 20;              // maximum number of older blocks for trace
int nbhist=-1;                      // current block for trace
uint64_t poptrace1[N2*nNhist];      // trace of first N*nNhist time points of population of N most frequent genes
uint64_t acttrace1[N2*nNhist];      // trace of first N*nNhist time points of activity of N most frequent genes
uint64_t acttraceq1[N2*nNhist];     // trace of first N*nNhist time points of activity of N most frequent quad patterns
unsigned char acttraceqt1[N2*nNhist];// type of entry in acttraceq : 1 quad, 0 smallpatt (<65536) i.e. corresponding to isnode
uint64_t working[N2];               // working space array for calculating genealogies and doing neighbour bit packing
int npopulation[N];                 // number of live sites (gol 1s) in last N time steps up to current population
int npopulation1[N*nNhist];         // number of live sites (gol 1s) in first N*nNhist time steps
int nspeciesgene,nallspecies;       // number of gene species in current population, and that have ever existed
int nallclones;                     // number of clones stored in hash table
int nspeciesquad,nallspeciesquad;   // number of quad species in current population, and that have ever existed
int nspeciessmall,nallspeciessmall; // number of small pattern species now, and that have ever existed
int histcumlogpattsize[log2N+1];    // histogram of patterns binned on log scale according to power of two side enclosing square
int histcumpixelssqrt[N+1];         // histogram of patterns binned on an integer sqrt scale according to number of pixels
//------------------------------------------------ arrays for connected component labelling and tracking ------------------------------------------------
const int NLM = N2>>2;              // maximum number of discrete components possible N2/4
short unsigned int label[N2];       // labels for pixels in connected component labelling
short unsigned int oldlabel[N2];    // previous time step labels for connected component labelling
short unsigned int labelcc[N2];     // label array to reassemble and display individual components
typedef struct equivrec {           // equivalence record
  short unsigned int pt;            // parent of record : integer index into eqv array at lower integer location
  short unsigned int rank;          // rank of equivalence class tree, may be used to accelerate union-find
  short unsigned int size;          // number of cells (lattice sites) with this label in current array
  short unsigned int overlaps;      // reserved for future use, filling record overall size to single long int (64bit)
} equivrec;
equivrec eqv[NLM];                  // equivalences between labels
typedef struct component {          // data structure for identified connected components in gol array
    short unsigned int N,S,W,E;     // rectangular bounds of rectangle containing the component
    short unsigned int lastrc;      // last occupied row/col used to trace components wrapping across the N-1 to 0 border
    short unsigned int label,log2n; // label index of component, enclosing square size n=2^log2n
    short unsigned int patt;        // for small components of size 4x4 pixels or less, the image is encoded directly in the pattern patt
    uint64_t quad;                  // hashkey for quadtree node of subimage for component
    unsigned int pixels;            // number of pixels on
    float gcolor  ;                 // inherited drifting color hue mixed from connected comp's (structure needs whole nr of 64-bit words for ndarray python comm.)
} component;
component complist[NLM];            // current array of components
component oldcomplist[NLM];         // old (previous time step) list of components
int ncomponents = 0;                // current number of components (= current number of labels)
int oldncomponents = 0;             // number of components in previous time step
typedef struct connection {         // connection of connected component at time t to another component at t-1
    short unsigned int oldlab;      // label at time t-1
    short unsigned int newlab;      // label at time t
    unsigned int next;              // next connection index in list of backward connections for a given labelled component at time t (newlab)
    unsigned int nextf;             // next connection index in list of forward connections for a given labelled component at time t-1 (oldlab)
    unsigned int overlap;           // sum of neighboring pixels (9-nbhd) between old and new component
    float woverlap;                 // weighted overlap out of all forward connections from old component i.e. overlap/sum(overlaps))
    float aoverlap;                 // weighted woverlap out of all backward connections from new component i.e. woverlap/sum(woverlaps))
} connection;
connection connections[N2];         // open memory reservoir of connection nodes to use in connection lists
unsigned int connlists[NLM];        // entry points for connection list corresponding to each labelled component
unsigned int connlistsf[NLM];       // entry points for forward connection list corresponding to each old labelled component
unsigned int connlen[NLM];          // lengths of connection lists for connected components to previous time ie from t to t-1
unsigned int connlenf[NLM];         // lengths of connection lists for connected components to current time i.e. from t-1 to t
unsigned int connpref[NLM];         // preferred backward connected component at time t-1 for each component at time t
unsigned int connpreff[NLM];        // preferred forward connected component at time t for each component at time t-1
int connused = 0;                   // used connection nodes
int gcolors = 0;                    // use inherited colors to color connected components
//............................................... optimal linear assignment t-1 to t ....................................................................
// #include "lapjv.h"               // modified from Tomas Kazmar python interfaced implementation of Jonker-Volgenant LAPMOD algorithm, if using lapmod.c
unsigned int iilap[NLM];            // indices of start of each variable length row in sparse cost matrix : first entry 0, last entry nclap
unsigned int cclap[N2];             // sparse cost matrix for mapping connected components at t-1 to t (overlaps) containing nclap entries
unsigned int recolor[NLM];          // recolor connected component based on colors of connected components and random drift
short unsigned int kklap[N2];       // column indices for successive entries in cost matrix
short unsigned int xlap[NLM];       // returned list of assignments: columns assigned to rows
short unsigned int ylap[NLM];       // returned list of assignments: rows assigned to columns
short unsigned int dist[NLM];       // distance along augmented paths for Hopcroft Karp matching algorithm: maxmatch
short unsigned int relabel[NLM];    // array to relabel connected components matching to be compatible with previous step
short unsigned int oldrelabel[NLM]; // old relabel array at previous time step
short unsigned int queue_array[NLM];// array for queue used in Hopcroft Karp matching algorithm: maxmatch
int  nlap;                          // number of connected components at t-1 entering into the assignment, i.e. n for LAPMOD
int  nclap;                         // number of edges between connected comp's t-1 to t entering into the assignment, == no. of cost matrix entries
int  nmatched;                      // number of matched old labels in current label set
//------------------------------------------------ planes and configuration offsets----------------------------------------------------------------------
int offdx=0,offdy=0,offdt=0;        // display chosen offsets for glider analysis with colorfunction 8
int Noff = 9;                       // number of offsets
int **offsets;                      // array of offsets (2D + time) for planes
int *histo;
int numHisto;
// initialize planes:
#define maxPlane 8                  /* maximum number of planes allowed : values 2,4,8 allowed */
int curPlane = 0;                   // current plane index
int newPlane = 1;                   // new plane index
int numPlane = maxPlane;            // number of planes must be power of 2 to allow efficient modulo plane
uint64_t *planes[maxPlane];         // ring buffer planes of gol array states
uint64_t *planesg[maxPlane];        // ring buffer planes of golg genes
uint64_t *planesgs[maxPlane];       // ring buffer planes of golgstatus bits
uint64_t *planesb[maxPlane];        // ring buffer planes of birth id for clone ids
uint64_t *planesr[maxPlane];        // ring buffer planes of birth id for resources
uint64_t plane0[N2];                // gol   0
uint64_t plane1[N2];                // gol   1
uint64_t planeg0[N2];               // golg  0
uint64_t planeg1[N2];               // golg  1
uint64_t planegs0[N2];              // golgs 0
uint64_t planegs1[N2];              // golgs 1
uint64_t planeb0[N2];               // golb  0
uint64_t planeb1[N2];               // golb  1
uint64_t planer0[N2];               // golr  0
uint64_t planer1[N2];               // golr  1
#if maxPlane > 2
uint64_t plane2[N2];                // gol   2
uint64_t plane3[N2];                // gol   3
uint64_t planeg2[N2];               // golg  2
uint64_t planeg3[N2];               // golg  3
uint64_t planegs2[N2];              // golgs 2
uint64_t planegs3[N2];              // golgs 3
uint64_t planeb2[N2];               // golb  2
uint64_t planeb3[N2];               // golb  3
uint64_t planer2[N2];               // golr  2
uint64_t planer3[N2];               // golr  3
#endif
#if maxPlane > 4
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
uint64_t planeb4[N2];               // golb  4
uint64_t planeb5[N2];               // golb  5
uint64_t planeb6[N2];               // golb  6
uint64_t planeb7[N2];               // golb  7
uint64_t planer4[N2];               // golr  4
uint64_t planer5[N2];               // golr  5
uint64_t planer6[N2];               // golr  6
uint64_t planer7[N2];               // golr  7
#endif
//------------------------------------------------------- fast macros for pattern counting and random number generator ---------------------------------
                                    // Wikipedia "Xorshift" rewritten here as inline macro &
                                    // Vigna, Sebastiano. "xorshift*/xorshift+ generators and the PRNG shootout". Retrieved 2014-10-25.
uint64_t randstate[2];              // State for xorshift pseudorandom number generation. The state must be seeded so that it is not zero
#define RAND128P(val) {                                                       \
    uint64_t x = randstate[0]; uint64_t y = randstate[1];                       \
    randstate[0] = y;    x ^= x << 23;  randstate[1] = x ^ y ^ (x >> 17) ^ (y >> 26);  \
    val = randstate[1] + y;}
int ranseed = 1234;
//.......................................................................................................................................................
const uint64_t m1  = 0x5555555555555555; //binary: 0101...           Constants for Hamming distance macro POPCOUNT64C
const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...
const uint64_t m10 = 0x1111111111111111; //binary: 00010001...       Additional constant for 4-bit byte distance macro POP4COUNT64C
#define POPCOUNT64C(x, val) {                  /* Wikipedia "Hamming Weight" popcount4c alg */  \
    uint64_t xxxx;                             /* define copy of x argument so that we do not change it */ \
    xxxx = x;                                  /* copy x argument */ \
    xxxx -= (xxxx >> 1) & m1;                  /* put count of each 2 bits into those 2 bits */ \
    xxxx = (xxxx & m2) + ((xxxx >> 2) & m2);   /* put count of each 4 bits into those 4 bits */ \
    xxxx = (xxxx + (xxxx >> 4)) & m4;          /* put count of each 8 bits into those 8 bits */ \
    val = (xxxx * h01) >> 56;}                 /* left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... */
#define POP4COUNT64C(x, val) {                 /* Count number of 4-bit bytes which are non zero */  \
    uint64_t xxxx;                             /* define copy of x argument so that we do not change it */ \
    xxxx = x;                                  /* copy x argument */ \
    xxxx = (xxxx | (xxxx >> 1)) & m1;          /* put or of each 2 bits into right bit of pair */ \
    xxxx = (xxxx | (xxxx >> 2)) & m10;         /* put or of each quartet bits into right bit of quartet */ \
    xxxx = (xxxx + (xxxx >> 4)) & m4;          /* put count of each 8 bits into those 8 bits */ \
    val = (xxxx * h01) >> 56;}                 /* left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... */
//.......................................................................................................................................................
#define PATTERN2_32(x, pat, found) {           /* find number of 2-bit aligned 2-bit copies of pattern pat in 32-bits of 64-bit integer x */ \
    const uint64_t r1_32 = 0x11111111ull;       \
    const uint64_t r3_32 = 0x33333333ull;       \
    const uint64_t r5_32 = 0x55555555ull;       \
    const uint64_t rf_32 = 0xffffffffull;       \
    uint64_t xxxx;                             /* define internal variable */ \
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
    uint64_t xxxx;                             /* define internal variable */ \
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
    uint64_t xxxx;                             /* define internal variable */ \
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
    uint64_t mmmm,mmmq;                        /* calculates position of rightmost (lsb) 1 : arguments must be of types uint64_t and int respectivley */ \
    int cccc;                                  /* takes on successive integer values 32,16,84,2,1 */ \
    int tttt;                                  /* logical to integer variable true=one false=zero : if a 1 in v under mask mmmm */ \
    c=v?0:1;                                   /* c will contain count of number of zeros on right of last one, here if v is all zeros then start from 1 */ \
    mmmm=~0ull;                                /* initially all ones, this is the mask from previous stage in loop below */ \
    for (cccc=1<<5;cccc>0;cccc>>=1) {          /* loop over cccc goes 32,16,8,4,2,1 : the amount of shift used in mask construction */ \
        mmmq = mmmm;                           /* query mask mmmq is to be the mask used to query if a one is under it at this stage, start with old mask */ \
        mmmq &= mmmm^(mmmm<<cccc);             /* mmmq: 0xffffffff00000000, 0xffff0000ffff0000, 0xff00ff00ff00ff00, 0xf0f0f0f0f0f0f0f0, 0xccc..., 0xaaa...*/ \
        tttt = v&mmmq?0:1;                     /* tttt is zero if a one under the query mask, one otherwise */ \
        mmmm=mmmq^(tttt*mmmm);                 /* the new mask for next stage is the query mask if a one is under it, otherwise the other half of mmmm */ \
        c+=tttt*cccc;                          /* the right zero counter is incremented by the length of the current interval cccc if a one was not under mask */ \
    }                                          /* note that Anderson's algorithm was incorrect, see also profile comparison in standalone lsb64.c */ \
}                                              /* this macro calculates the LSB 1 (ie from the bottom) not the MSB 1 (ie from the top) that the integer log function finds. */
//.......................................................................................................................................................
#define membrane (((ij>>log2N)==((N>>1)-(initfield>>1)-1) || ((ij>>log2N)==((N>>1)-(initfield>>1)-2))) && (ij & 0x1) ? 2 : 1) /* formula for membrane of death */
//----------------------------------------------------- list of subroutines -----------------------------------------------------------------------------
//......................................................  fast integer processing   .....................................................................
// integerSqrt          direct bit processing algorithm to implement integer sqrt (largest integer smaller than sqrt) : but floating point sqrt is faster
// log2r                fast integer logarithm working only for arguments which are powers of 2 (not used but slightly faster than log2a when applicable)
// log2lower            fast integer logarithm working for all integers : largest integer smaller than or equal to logarithm base 2 of argument
// log2upper            fast integer logarithm working for all integers : smallest integer larger than or equal to logarithm base 2 of argument
// sqrtupper            fast integer sqrt working for all integers : smallest integer larger than or equal to sqrt of argument
// randprob             random event with probability determined by a 32 bit unsigned integer iprob as iprob / 2^32 using RAND128 and uint64_t
//......................................................  color functions   .............................................................................
// set_color            inline function to assign rainbow colors
// label_color          map label (quasi-uniformly) into larger unsigned colour space
// rgba                 converts hue (0..1) to rgb+alpha
// mix_color            mix colors from overlapping components, add random drift of colour
// delay                time delay in ms for graphics
// printxy              terminal screen print of array on xterm
// colorgenes           colour display of genes in one of 12 colorfunction modes, including activities, pattern analysis, genealogies and glider detection
//......................................................  selection of genes for birth  .................................................................
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
//...................................................... hash table management for genes and clones .....................................................
// hashaddgene          add new gene to hash table, increment popln and ancestor information for genes already encountered. If mutation calls hashaddclone
// hashdeletegene       decrements population count for gene, printing error message if not found or already zero. Calls hashdeletefromclone
// hashgeneextinction   record gene extinctions in hash gene table, counting number of extinctions
// hashgeneactivity     update activity of gene
// hashaddclone         add new clone to hash table, increment popln information for existing clones
// hashdeletefromclone  decrements population count for clone, printing error message if not found or already zero
// hashcloneactivity    update activity of clone
//......................................................  pattern storage and analysis  .................................................................
// patt_hash            hash function for a pattern specified by 4 64-bit (8x8) patterns
// newkey               assign one of free keys kept in a linked list, allocated 1024 at a time (used if hash-key already occupied with different quadtree)
// node_hash            hash function for a node specified by 4 64-bit pointers
// hash_patt16_store    store new pattern in small pattern table
// hash_patt16_find     find quadtree hash for pattern (leaf of quadtree consists of 4 64bit integers defining a 16x16 pixel array)
// hash_node_store      store new node with hashkey h and subnodes nw,ne,sw,se in hash table
// hash_node_find       find quadtree hash for node (node is specified by its four quadrant pointers (64 bit))
// quadimage            construct quadtree for an entire image or connected component, reporting if the image has been found previously, returning hashkey
// labelimage           rebuild image in a chosen label array from quadimage at chosen offset with chosen label
//......................................................  neighborhood processing  ......................................................................
// pack012neighbors     pack all up to 2nd neighbours in single word
// pack0123neighbors    pack all up to 3rd neighbours in single word
// pack49neighbors      fast routine to pack all up to 1st,2nd,3rd neighbours in single word : order of bits dictated by hierarchical assembly
// pack16neighbors      pack 4x4 blocks in single uint64_t word using 16 bits
// unpack16neighbors    unpack 16 bit word to 4x4 block at offset in full array of labels, marking with chosen label
// log2size             return log2 of linear size of pattern in integer power of 2 for small patterns 0-65535
// pack64neighbors      pack 8x8 blocks in single uint64_t (long) word with 64 bits
// unpack64neighbors    unpack 64 bit word to 8x8 block at offset in full array of labels, marking with chosen label
// compare_neighbors    compare packed pack neighbours with one given x,y shift of arbitrary size
// compare_all_neighbors compare packed pack neighbours with all nearest neighbour x,y shifts
// packandcompare       pack and compare either all 1-shifted 3-neighbourhoods with t=-1 or chosen (dx,dy,dt) 3-neighbourhoods
//......................................................  fast component labelling  ......................................................................
// lab_union            disjoint rank union of equivalence classes returning common root
// label_cell           label a cell (site) in the cellular automata with first pass label of connected component
// label_cell_Wu        alternative version according to secondary paper by Wu discussing Suzuki decision tree
// checklabels          check that the label tree consistently points to labels of lower values as we go via parents to root
// flattenlabels        flatten label tree so that each label points to its unique root
// label_components     do two-pass fast component labelling with 8-neighbour using Suzuki decision tree, rank union and periodic BCs, connect t-1 labels with t
// extract_components   extract labelled components to list of subimages embedded in square of side 2^n, each stored in a quadtree hash table
//........................................................  simulation update for different symmetries  ..................................................
// update_23            update gol, golg, golgstats for a single synchronous time step : for selection 0-7 with fixed GoL rule departures in repscheme
// update_lut_sum       update version for gene encoding look up table for totalistic survival and birth (disallowing 0 live neighbour entries) sel 8,9
// update_lut_dist      update version for gene encoding look up table for survival and birth based on corner & edge sums (2*19 states disallowing s=0,1,7,8): sel 10,11
// update_lut_canon_rot update version for gene encoding look up table for canonical rotation survival and birth (2*32 states, disallowing 0,1,7,8 entries) : sel 12,13
// update_lut_2D_sym    update version all different configurations under the standard 2D 4-rotation and 4-reflection symmetries are distinguished: sel 14,15
// genelife_update      master routine to call specific model update, collect statistics if required and rotate planes
//........................................................  initialization and file IO  ...................................................................
// initialize_planes    initialize periodic sequence of planes to record rolling time window of up to maxPlanes time points (≤8)
// readFile             read file of gol/golg array (32x32) data
// writeFile            write file of gol/golg array (32x32) data
// testmacros           test macros used to accelerate processing (usually not called): FIRST1INDEX, PATTERN4, PATTERN8
// initialize           initialize simulation parameters and arrays
//.........................................................  set from python driver  ......................................................................
// set_colorfunction    set color function integer from GUI for use in patterning and coloring display
// setget_act_ymax      set activity ymax for scaling of gene activity plot
// setget_act_ymaxq     set activity ymax for scaling of quad activity plot
// set_selectedgene     set selected gene for highlighting from current mouse selection in graphics window
// set_offsets          set offsets for detection of glider structures in display for color function 8
// set_quadrant         set the pair of bits in repscheme (or survivalmask or overwritemask) used for quadrant variation 0-6
// set_randominflux     change the randominflux activation for rbackground if nonzero or continual updating of init field with random states and genes : 2,1,0
// set_rbackground      set the backround random live gene input rate per frame and site to rbackground/32768
// set_repscheme_bits   set the two of the repscheme (or survivalmask or overwritemask) bits corresponding to the selected quadrant
// set_repscheme        set repscheme from python
// set_rulemod          set rulemod from python
// set_surviveover      set the two masks for survival and overwrite from python (survivalmask, overwritemask)
// set_vscrolling       set vertical scrolling to track fronts of growth in vertical upwards direction
// set_noveltyfilter    set novelty filter for darkening already encountered components in connected component display (colorfunction 9)
// set_activity_size_colormode set colormode by size for colorfunction 10 : 0 by ID  1 log2 enclosing square size 2 use #pixels 3 use sqrt(#pixels)
// set_gcolors          set connected component colors as inherited colors from colliding connected components with random drift
// set_seed             set random number seed
// set_nbhist           set nbhist N-block of time points for trace from GUI for use in activity and population display traces
// set_genealogycoldepth set genealogycoldepth for colorfunction=11 display
// set_ancestortype     set ancestortype for genealogy display and return of first (0), clonal (1) or first & clonal in 2 windows (2)
// set_stash            stash current gol,golg in stashgol, stshgolg
// set_info_transfer_h  set information transfer histogram display value (0,1) from python
// set_activityfnlut    set collection of functional activity statistics corresponding to functional aggregate of genes by non-neutral bits
// set_colorupdate1     control update of colorgenes and regular print statements via flag colorupdate1
// set_colorfunction2   choice of colorfunction for window 2
//..........................................................  get to python driver  .....................................................................
// get_stash            retrieve current gol,golg from stashed values
// get_log2N            get the current log2N value from C to python
// get_curgol           get current gol array from C to python
// get_curgolg          get current golg array from C to python
// get_curgolbr         get current golb and golr arrays from C to python
// get_stats            get the traced statistics from C to python
// get_acttrace         get current acttrace array C to python
// get_poptrace         get current poptrace array C to python
// get genealogytrace   get current trace of genealogies to python
// get_nspecies         get number of species from C to python
// get_genealogydepth   get depth of genealogies returned by get_genealogies()
// get_genealogies      // get current population genealogies (currently near end of code)
// get_hist             get the histogram from C to python
// get_activities       get the current activity statistics of genes from C to python
// get_all_activities   get all activity statistics of genes (since t=0) from C to python
// get_quad_activities  get current activity statistics of quads (since t=0) from C to python
// get_small_activities  get current activity statistics of smallpats (since t=0) from C to python
// get_all_quad_activities  get all activity statistics of quads (since t=0) from C to python
// get_all_small_activities  get all activity statistics of smallpats (since t=0) from C to python
// get_components       get all current connected component data structures
// get_smallpatts       get array of small pattern data structures including sizes and activities
// get_quadnodes        get all hashed quadnodes including hashkey, sizes and activities
// get_genes            get all hashed genes with data structures including activity counts, extinctions etc
// get_curgolgstats     get current golgstats array C to python
// get_sorted_popln_act return sorted population and activities (sorted by current population numbers)
//..........................................................  comparison functions  ....................................................................
// cmpfunc              compare gene values as numerical unsigned numbers
// cmpfunc1             compare gene counts in population
// cmpfunc2             compare gene values corresponding to given number index in hash table
// cmpfunc3             compare population counts of hash stored genes
// cmpfunc3c            compare population counts of hash stored clones
// cmpfunc3q            compare pixel counts (pop1s) of hash stored quad patterns
// cmpfunc3qs           compare pixel counts (pop1s) of hash stored small patterns
// cmpfunc4             compare birth times of hash stored genes
// cmpfunc5             compare common genealogy level gene values of hash stored genes
// cmpfunc5c            compare common clonealogy level clone values of hash stored clones
// cmpfunct6            compare according to ancestry in genealogytrace using activity ordering
// cmpfunc7             compare according to ancestry in genealogytrace using population size ordering
//..........................................................  gene and pattern analysis of dynamics .....................................................
// countconfigs         count the configs with python specified offsets in (x,y,t)
// tracestats           record the current stats in time trace
// countspecies1        count genes with gene array specified as input parameters
// countspecies         count different genes with genes specified at current time point
// countspecieshash     count different genes in current population from record of all species that have existed
// totalpoptrace        calculates and returns current population size and store in scrolling population trace array npopulation
// genefnindex          calculate index based on bits masked in from survival and birth masks only
// activitieshash       calculate array of current gene activities and update acttrace array of genes in activity plot format
// activitieshashquad   calculate array of current quad activities and update acttraceq array of patterns in activity plot format
// get_genealogies      calculate and retrieve or display genealogies, depending on size of array passed (0: display only, >0: retrieve only)
// clonealogies         calculate and display clonealogies: genealogies of clones
// get_gliderinfo       get information about gliders from packed array representation
//-------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------- mathematical fns ----------------------------------------------------------------------
extern inline int integerSqrt(int n) {                  // the largest integer smaller than the square root of n (n>=0)
    int shift,nShifted,result,candidateResult;
    // only works for n >= 0;
    // Find greatest shift.
    shift = 2;
    nShifted = n >> shift;

    while (nShifted != 0) {
        shift += 2;
        nShifted >>= 2;
    }
    shift -= 2;

    // Find digits of result.
    result = 0;
    while (shift >= 0) {
        result <<=  1;
        candidateResult = result + 1;
        if (candidateResult*candidateResult <= (n >> shift))
            result = candidateResult;
        shift -= 2;
    }
    return result;
}
//.......................................................................................................................................................
extern inline unsigned int log2r(unsigned int v) { // find the log2 of v = power of 2, warning wrong answers for v not power of 2
    // From:  https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
    static const unsigned int b[] = {0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0,
                                     0xFF00FF00, 0xFFFF0000};
    register unsigned int r = (v & b[0]) != 0;
    r |= ((v & b[4]) != 0) << 4;
    r |= ((v & b[3]) != 0) << 3;
    r |= ((v & b[2]) != 0) << 2;
    r |= ((v & b[1]) != 0) << 1;
    return(r);
}
//.......................................................................................................................................................
extern inline unsigned int log2lower(unsigned int v) { // find the integer log2 of v (lower): works for all 32 bit integer v > 0
// Adapted by John McCaskill starting from the power of two version at https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
// log2upper uses the fact that this routine gives log2lower(0)=0

    static const unsigned int b[] = {0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0,
                                     0xFF00FF00, 0xFFFF0000};
    register int k = 4;
    register unsigned int r,s;
    if(!v) return 0;
    s = ((v & b[k])   != 0); r  = s; v &= s ? b[k] : ~b[k]; r<<=1;
    s = ((v & b[--k]) != 0); r |= s; v &= s ? b[k] : ~b[k]; r<<=1;
    s = ((v & b[--k]) != 0); r |= s; v &= s ? b[k] : ~b[k]; r<<=1;
    s = ((v & b[--k]) != 0); r |= s; v &= s ? b[k] : ~b[k]; r<<=1;
    s = ((v & b[--k]) != 0); r |= s;

    return(r);
}
//.......................................................................................................................................................
extern inline unsigned int log2upper(unsigned int v) { // find the integer log2 of v (upper): works for all 32 bit integer v > 0

    return v ? 1 + log2lower(v-1) : 0;              // assuming as in this implementation that log2lower(0) == 0
}
//.......................................................................................................................................................
extern inline unsigned int sqrtupper(unsigned int v) { // smallest integer larger than or equal to sqrt of argument

    return v ? 1 + (int) sqrt(v-1) : 0;
}
//..............................................................  randprob  .............................................................................
extern inline uint64_t randprob(unsigned int uprob, unsigned int randnr) {
    return(randnr < uprob ? 1ull : 0ull);
}
//--------------------------------------------------------------- colorgenes ----------------------------------------------------------------------------
extern inline void setcolor(unsigned int *color,int n) { // for coloring by quad size...
    // rainbow colors from running from fields::tim.colors(11) in R by Tim Hoar : assuming log2N <= 10 (limited by display size)
    // "#00008F" "#0000F5" "#005AFF" "#00BDFF" "#23FFDC" "#87FF78" "#ECFF13" "#FFAD00" "#FF4A00" "#E40000" "#800000"
    static const unsigned int rainbow[]={0x800000FF, 0xE40000FF, 0xFF4A00FF, 0xFFAD00FF, 0xECFF13FF, 0x87FF78FF, 0x23FFDCFF, 0x00BDFFFF, 0x005AFFFF, 0x0000F5FF, 0x00008FFF};

    unsigned int mycol;
    mycol = rainbow[n];
    color[2] = (mycol >> 24) & 0xff; // B
    color[1] = (mycol >> 16) & 0xff; // G
    color[0] = (mycol >> 8) & 0xff;  // R
}
//.......................................................................................................................................................
extern inline unsigned int labelcolor( short unsigned int label) {
    uint64_t mask;
    unsigned int color;
    mask = (uint64_t) label;
    mask = mask * 11400714819323198549ull;          // map label (quasi-uniformly) into larger unsigned colour space
    color = mask >> (64 - 32);                      // hash with optimal prime multiplicator down to 32 bits
    color |= 0x080808ffull;                         // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
    return color;
}
//.......................................................................................................................................................
extern inline unsigned int rgba( float hue) {                            // converts hue (0..1) to rgb+alpha
    unsigned int r,g,b;
    unsigned int k;
    float hue6 = hue*6.0;

    k = (unsigned int) (6.*hue); if (k==6) k = 0;  // deal with possible round off error
    r=g=b=0x00;
    switch (k) {
        case 0: r = 0xff; g = (unsigned int)((hue6)*256.0);   g = g>0xff?0xff:g;break;
        case 1: g = 0xff; r = (unsigned int)((2.-hue6)*256.0);r = r>0xff?0xff:r;break;
        case 2: g = 0xff; b = (unsigned int)((hue6-2.)*256.0);b = b>0xff?0xff:b;break;
        case 3: b = 0xff; g = (unsigned int)((4.-hue6)*256.0);g = g>0xff?0xff:g;break;
        case 4: b = 0xff; r = (unsigned int)((hue6-4.)*256.0);r = r>0xff?0xff:r;break;
        case 5: r = 0xff; b = (unsigned int)((6.-hue6)*256.0);b = b>0xff?0xff:b;break;
        default: fprintf(stderr,"step %d Error in hut to rgb color switch - case %d out of bounds\n",totsteps,k);
    }
    return ((r<<8) | (g<<16) | (b<<24) | 0xff);
}
//.......................................................................................................................................................
extern inline float mixcolor( short unsigned int label,uint64_t rand) { // mix colors from overlapping components, add random drift of colour
    unsigned int conn;
    float color,color1,x1,y1,x,y;
    float eps = 0.0001*gcolors;
    #define PI 3.14159265
    const float i2pi = 1./(2.*PI);

    color = x = y = 0.0;
    conn = connlists[label];                                            // NB conn is not a label but an index in the array of possible connections
    while(conn) {
        color1= oldcomplist[connections[conn].oldlab].gcolor;           // mix colors based on aoverlap weights
        x1 = cosf(2.*PI*color1);                                        // averaging of circular variable requires converting to x,y vector and average
        y1 = sinf(2.*PI*color1);
        /* color2 = atan2f(y1,x1)*i2pi;
        color2 = color2 < 0. ? color2 + 1. : color2;
        if (fabsf(color2-color1) > eps) fprintf(stderr,"step %d label %d non inverted arctan color1 %f color2 %f\n",totsteps,label,color1,color2);*/
        x += x1*connections[conn].aoverlap;
        y += y1*connections[conn].aoverlap;
        conn=connections[conn].next;
    }
    color = atan2f(y,x)*i2pi;                                           // atan2 returns principal value in range [-PI, PI]
    color = color < 0. ? color + 1. : color;                            // for negative values add 2*PI to get back to range [0,2PI]
    color += eps*(((float)(rand&0xff)) - 127.5);                        // random drift of color
    color = color>1.0 ? color-1.0 : color;
    color = color<0.0 ? color+1.0 : color;
    return color;
}
//.......................................................................................................................................................
void delay(int milliseconds) {
    long pause;
    clock_t now,then;
    pause = milliseconds*(CLOCKS_PER_SEC/1000);
    now = then = clock();
    while( (now-then) < pause )
        now = clock();
}
//.......................................................................................................................................................
void printxy (uint64_t gol[],uint64_t golg[]) {                         // print the game of life configuration
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
//.......................................................................................................................................................
extern inline void golr_digest (uint64_t golr, unsigned int *mismatchmin, unsigned int *mismatchmax, unsigned int *period, int *pershx, int *pershy) {
                                                        // extract optimal period match: min, max numbers of mismatches,period and x,y offsets for opt. match period
    uint64_t gdiff;
    int j,k,d,d0,d1,dx,dy;
    unsigned int jper;
    int nbx[8] = {-1,0,1,1,1,0,-1,-1};
    int nby[8] = {-1,-1,-1,0,1,1,1,0};
    int dsx,dsy;

    
    gdiff = golr;                                       // variable golr holds dynamical record
    d0 = 16; d1 = 0;                                    // min,max number of mismatches

    jper = 0;
    for (j=0;j<15;j++) {
        gdiff = (gdiff>>4)|((gdiff&0xfull)<<60);        // rotate record cyclically by one time step
        POP4COUNT64C((golr^gdiff),d);                    // number of difference positions between gene and gdiff
        if(d<d0) {
            d0=d;
            jper = j;
        }
        if(d>d1) d1=d;
    }
    *mismatchmin = d0;
    *mismatchmax = d1;
    *period = jper;
    
    dsx = dsy = dx = dy = 0;
    for (j=0;j<16;j++) {
        gdiff = (golr>>(j<<2))&0xfull;
        if(gdiff&0x8ull) {
            k = ((unsigned int) gdiff) & 0x7;
            dx += nbx[k];
            dy += nby[k];
        }
        if (j>=jper) {
            dsx+=dx;
            dsy+=dy;
            dx = dy = 0;
        }
    }
    // average periodic displacements are now dsx/(16-jper), dsy/(16-jper)
    *pershx = dsx;
    *pershy = dsy;
}
//.......................................................................................................................................................
void colorgenes( int cgolg[], int NN2, int colorfunction, int winnr, int nfrstep) {
    uint64_t gene, gdiff, g2c, mask, quad, clone;
    uint64_t *ggol,*ggolg,*ggolb,*ggolr,*ggolgstats;
    int ij,k,j,nfrPlane,nbeven,activity,popcount,labelxy;
    unsigned int d,d0,d1,d2,jper;
    unsigned int color[3],colormax;
    double rescalecolor;
    uint64_t * traceptr;
    int *npopptr;
    quadnode *q;
    static int numones[16]={0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};
    int labelimage(uint64_t hashkeypatt, short unsigned int labelimg[], short unsigned int label, int offset);
    extern inline int log2size(const short unsigned int golpw);
    void get_gliderinfo(uint64_t outgliderinfo[], int narraysize);
    static int firstdebug = 1;
    
    if (nfrstep==0) {
        ggol  = gol;
        ggolg = golg;
        ggolb = golb;
        ggolr = golr;
        ggolgstats = golgstats;
    }
    else {
        nfrPlane = (curPlane+nfrstep) % numPlane;
        ggol = planes[nfrPlane];                                               // get planes of gol,golg,golb,golr,golgstats data
        ggolg = planesg[nfrPlane];
        ggolgstats = planesgs[nfrPlane];
        ggolb = planesb[nfrPlane];
        ggolr = planesr[nfrPlane];
    }

    if(colorfunction==0) { // colorfunction based on multiplicative hash
        // see https://stackoverflow.com/questions/6943493/hash-table-with-64-bit-values-as-key/33871291
        for (ij=0; ij<N2; ij++) {
            if (ggol[ij]) {
                gene = ggolg[ij];
                if (gene == 0ull) gene = 11778ull; // random color for gene==0
                // mask = (gene * 11400714819323198549ul) >> (64 - 8);   // hash with optimal prime multiplicator down to 8 bits
                // mask = (gene * 11400714819323198549ul) >> (64 - 32);  // hash with optimal prime multiplicator down to 32 bits
                mask = gene * 11400714819323198549ull;
                mask = mask >> (64 - 32);   // hash with optimal prime multiplicator down to 32 bits
                if(rulemod & 0x4) mask = (((ij>>log2N)==((N>>1)-(initfield>>1)-1) || ((ij>>log2N)==((N>>1)-(initfield>>1)-2))) && (ij & 0x1) ? 0x00ffffffull : mask);
                mask |= 0x080808ffull; // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                cgolg[ij] = (int) mask;
            }
            else {
                if(rulemod & 0x4) mask = (((ij>>log2N)==((N>>1)-(initfield>>1)-1) || ((ij>>log2N)==((N>>1)-(initfield>>1)-2))) && (ij & 0x1) ? 0x00ffffffull : 0ull);
                else mask = 0ull;
                cgolg[ij] = (int) mask;
            }
        }
    }
    else if(colorfunction<4){
        for (ij=0; ij<N2; ij++) {
            if (gol[ij]) {
                gene = ggolg[ij];
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
                        case 8 :                                            // colors according to nr of non zero LUT entries d from blue to red in increasing d, except for repselect 7
                                 if (((repscheme>>4) & 0x7) == 7) {d = d & 0x3; mask = d==3 ? 0xf0f0f0ff : ((0xff<<(d<<3))<<8)+0xff;} // repselect 7 : sc-st-we-pa: red-green-blue-white
                                 else if (ncoding==1) {gdiff=gene&0xffff;POPCOUNT64C(gdiff,d);mask = d==16 ? 0x0000ffff : ((((15-d)<<20)+(abs((int)d-8)<<12)+(d<<4))<<8) + 0xff;}
                                 else if (ncoding==2) {PATTERN2_32(gene,0x3,d);mask = d==16 ? 0x00ffffff : ((((15-d)<<20)+(abs((int)d-8)<<12)+(d<<4))<<8) + 0xff;}
                                 else {PATTERN4(gene,0xf,d);mask = d==16 ? 0x00ffffff : ((((15-d)<<20)+(abs((int)d-8)<<12)+(d<<4))<<8) + 0xff;}
                                 break;
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
                    if(ggolgstats[ij]&F_nongolchg) mask = 0x00ffffff;        // color states changed by non GoL rule yellow
                }
                else if (colorfunction==3) {
                    if(ggolgstats[ij]&F_survmut) mask = 0xff00ffff;  // color states for surviving fresh mutants (non-replicated) purple/pink
                    else if(ggolgstats[ij]&F_notgolrul) mask = 0x00ffffff;  // color states for not GoL rule yellow
                }

                cgolg[ij] = (int) mask;
            }
            else cgolg[ij] = 0;
        }
    }
    else if(colorfunction==4){                                  // activities
        const int colpopmax = 225;                                      // need to bring this parameter up to python
        const double colpoplow = 50.0;
        const double logcolpopmax= log(colpoplow+(double) colpopmax);

        if(nbhist==-1) {
            traceptr =&acttrace[0];
            npopptr = &npopulation[0];
        }
        else {
            traceptr = &acttrace1[N2*(nbhist>>1)];
            npopptr = &npopulation1[N*(nbhist>>1)];
        }
        nbeven=(nbhist==-1)?1:1-(nbhist&0x1);
        for (ij=0; ij<N2; ij++) {
            int i,i1;
            i=ij&Nmask;
            i1=nbeven?0:(i<(N>>1)?N>>1:N2-(N>>1));
            gene=traceptr[ij+i1];
            if (gene == rootgene) mask = 0x3f3f3fff;            // grey color for background, all root genes
            else {
                if (gene == 0ull) gene = 11778L;                // random color for gene==0
                mask = gene * 11400714819323198549ul;
                mask = mask >> (64 - 32);                       // hash with optimal prime multiplicator down to 32 bits
                mask |= 0x080808ffull;                          // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                if(colpopmax && (diagnostics & diag_hash_genes)) {
                    popcount=0;
                    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) popcount = genedataptr->popcount;
                    else fprintf(stderr,"gene not found in colorfunction for activities\n");
                    if(popcount>colpopmax) popcount=colpopmax;
                    colormax=0;
                    for(d=0;d<3;d++) if((color[d]=( (mask>>(8+(d<<3))) & 0xff))>colormax) colormax=color[d];
                    if(popcount)
                        rescalecolor=(log(colpoplow+(double)popcount)/logcolpopmax)*((double)0xff/(double)colormax);
                    else
                        rescalecolor=0.2;
                    for(d=0;d<3;d++) color[d]=(unsigned int) (((double) color[d])*rescalecolor);
                    for(d=0,mask=0xff;d<3;d++) mask |= color[d]<<((d<<3)+8);
                }
            }
            i1=nbeven?0:(i<(N>>1)?N>>1:N-(N>>1));
            i=nbeven?0:(i<(N>>1)?0:N);
            if ((npopptr[i+((ij+i1)&Nmask)]>>log2N)==(N-1-(ij>>log2N))) mask = 0xffffffff;  // overlay plot with trace of density in white (except if pop=N2)
            cgolg[ij]= (int) mask;
        }
    }
    else if(colorfunction==5){                                  // populations of genes or clones
        int actmax = 0;                                         // need to bring this parameter up to python
        if(nbhist==-1) {
            traceptr =&poptrace[0];
            npopptr = &npopulation[0];
        }
        else {
            traceptr = &poptrace1[N2*(nbhist>>1)];
            npopptr = &npopulation1[N*(nbhist>>1)];
        }
        nbeven=(nbhist==-1)?1:1-(nbhist&0x1);
        for (ij=0; ij<N2; ij++) {
            int i,i1;
            i=ij&Nmask;
            i1=nbeven?0:(i<(N>>1)?N>>1:N2-(N>>1));
            gene=traceptr[ij+i1];
            if (gene == rootgene) mask = 0x3f3f3fff;            // grey color for background, all root genes
            else {
                if (gene == 0ull) gene = 11778L;                // random color for gene==0
                mask = gene * 11400714819323198549ul;
                mask = mask >> (64 - 32);                       // hash with optimal prime multiplicator down to 32 bits
                mask |= 0x080808ffull;                          // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                if(actmax && (diagnostics & diag_hash_genes)) {
                    popcount=0;
                    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) popcount = genedataptr->popcount;
                    else fprintf(stderr,"gene not found in colorfunction for population\n");
                    if(popcount>actmax) popcount=actmax;
                    colormax=0;
                    for(d=0;d<3;d++) if((color[d]=( (mask>>(8+(d<<3))) & 0xff))>colormax) colormax=color[d];
                    rescalecolor=(log((double)popcount)/log((double)actmax))*((double)0xff/(double)colormax);
                    for(d=0;d<3;d++) color[d]=(unsigned int) (((double) color[d])*rescalecolor);
                    for(d=0,mask=0xff;d<3;d++) mask |= color[d]<<((d<<3)+8);
                }
            }
            i1=nbeven?0:(i<(N>>1)?N>>1:N-(N>>1));
            i=nbeven?0:(i<(N>>1)?0:N);
            if ((npopptr[i+((ij+i1)&Nmask)]>>log2N)==(N-1-(ij>>log2N))) mask = 0xffffffff;  // overlay plot with trace of density in white (except if pop=N2)
            cgolg[ij]= (int) mask;
        }
    }
    else if(colorfunction==6 || colorfunction == 7){                    //genealogies
        int ancestortypec;
        if(colorfunction==6) k=2; else k=0;                             // double row for each ancestral step to make display more readable in colorfunction 6 mode
        if((colorfunction2==-1) || (ancestortype !=2)) ancestortypec = ancestortype;
        else if (winnr) ancestortypec = 1;
        else ancestortypec = 0;
        for (ij=0; ij<N2; ij++) {
            int i=ij&Nmask; int j=ij>>(log2N+k);
            if((ancestortypec==1) && (diagnostics &diag_hash_clones)) {
                clone=clonealogytrace[i+(j<<log2N)];
                activity = 0;
                if (clone == rootclone)   mask = 0x000000ff;            // black color for root
                else {
                    if (clone == 0ull) clone = 11778L;                  // random color for clone==0
                    mask = clone * 11400714819323198549ul;
                    mask = mask >> (64 - 32);                           // hash with optimal prime multiplicator down to 32 bits
                    mask |= 0x080808ffull; // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                    if((colorfunction==7) && (diagnostics & diag_hash_clones)) {                         // rescale color brightness by activity/activitymax
                        activity = 0;
                        if((clonedataptr = (clonedata *) hashtable_find(&clonetable, clone)) != NULL) activity = clonedataptr->activity;
                        else fprintf(stderr,"clone not found in colorfunction for clonealogy\n");
                        colormax=0;
                        for(d=0;d<3;d++) if((color[d]=( (mask>>(8+(d<<3))) & 0xff))>colormax) colormax=color[d];
                        rescalecolor=0.25+0.75*((double)(activity*255))/((double)(activitymax*colormax));         // integer version doesn't work
                        for(d=0;d<3;d++) color[d]=(unsigned int) (((double) color[d])*rescalecolor);     // rescale colors by activity/activitymax
                        for(d=0,mask=0xff;d<3;d++) mask |= color[d]<<((d<<3)+8);
                    }
                    else {};    // no rescaling of colours for case 6 (formerly cases 5 and 6)
                }
            }
            else {
                gene=genealogytrace[i+(j<<log2N)];
                activity = 0;
                if (gene == selectedgene)    mask = 0xffffffff;         // white color for selected gene
                else if (gene == rootgene)   mask = 0x000000ff;         // black color for root
                else if (gene == generepeat) mask = 0x3f3f3fff;         // grey color for repeated gene
                else {
                    if (gene == 0ull) gene = 11778L;                    // random color for gene==0
                    mask = gene * 11400714819323198549ul;
                    mask = mask >> (64 - 32);                           // hash with optimal prime multiplicator down to 32 bits
                    mask |= 0x080808ffull; // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                    if((colorfunction==7) && (diagnostics & diag_hash_genes)) {                           // rescale color brightness by activity/activitymax
                        activity = 0;
                        if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) activity = genedataptr->activity;
                        else fprintf(stderr,"gene not found in colorfunction for genealogy\n");
                        colormax=0;
                        for(d=0;d<3;d++) if((color[d]=( (mask>>(8+(d<<3))) & 0xff))>colormax) colormax=color[d];
                        rescalecolor=0.25+0.75*((double)(activity*255))/((double)(activitymax*colormax));         // integer version doesn't work
                        for(d=0;d<3;d++) color[d]=(unsigned int) (((double) color[d])*rescalecolor);     // rescale colors by activity/activitymax
                        for(d=0,mask=0xff;d<3;d++) mask |= color[d]<<((d<<3)+8);
                    }
                    else {};    // no rescaling of colours for case 6 (formerly cases 5 and 6)
                }
            }
            cgolg[ij]=(int) mask;
        }
    }
    else if(colorfunction == 8) {         // colorfunction based on packed bit pattern with multiplicative hash, compare with offset
        int it_nbhood2;
        it_nbhood2 = it_nbhood*it_nbhood+2;
        for (ij=0; ij<N2; ij++) {
                gene = golmix[ij];
                if(offdx==0 && offdy==0 && offdt==0) {
                    for (mask=0,k=0;k<8;k++) {
                        d1 = (gene>>(k<<3))&0xff;                   // no of mismatches : 0:0 black space perfect match is coded as 0xff=255
                        d1 = (d1 >= it_nbhood2-1) ? 0 : (it_nbhood2-1-d1);               // no of matches : 0:0 match set to 0 (no match) here
                        d1 = (d1 < it_nbhood2-1-15) ? 0 : d1-(it_nbhood2-1-15);          // 0 to 15 : perfect match is 15  (4 bits)
                        d1 = (d1==0xf) ? 0x1f : d1;                 // perfect match separated to value 31 (5 bits) for better contrast
                        if(k<3) mask+=d1<<(3+(k<<3));               // perfect match has full intensity colour
                        else if (k==3 && d1==0x1f) mask = (d1<<3)+(d1<<11)+(d1<<19); // the fourth channel has white colour : no others shown
                        else if (k<7 && d1==0x1f) mask+= (d1<<(3+((k-4)<<3)))+(d1<<(3+((k<6?k-3:0)<<3))); // mixed colours for NE SE SW
                        else if(d1==0x1f) mask+= (d1<<3)+(d1<<10)+(d1<<18); // mixed colour for NW
                    }
                    mask = (mask<<8)+0xff;
                }
                else {
                    POPCOUNT64C(gene,d);                            // assigns number of ones in gene to d. This 3 line version for one offset comparison
                    d=(d==64)?0:63-d;
                    mask = (d==63) ? 0xffffffff : ((((d&3)<<22)+(((d>>2)&3)<<14)+(((d>>4)&3)<<6))<<8) + 0xff;
                }

                cgolg[ij] = (int) mask;
        }
        if(info_transfer_h) {                                       // display histograms of glider matching in eight directions N E S W NE SE SW NW
            
            uint64_t *binomialp;
            uint64_t binomial9[11] = {1, 9, 36, 84, 126, 126, 84, 36, 9, 1, 1};
            uint64_t binomial25[27] = {1, 25, 300, 2300, 12650, 53130, 177100, 480700, 1081575, 2042975, 3268760, 4457400, 5200300, 5200300, 4457400, 3268760, 2042975, 1081575, 480700, 177100, 53130, 12650, 2300, 300, 25, 1, 1};
            uint64_t binomial49[51] = {1, 49, 1176, 18424, 211876, 1906884, 13983816, 85900584, 450978066, 2054455634, 8217822536, 29135916264, 92263734836, 262596783764, 675248872536, 1575580702584, 3348108992991, 6499270398159, 11554258485616, 18851684897584, 28277527346376, 39049918716424, 49699896548176, 58343356817424, 63205303218876, 63205303218876, 58343356817424, 49699896548176, 39049918716424, 28277527346376, 18851684897584, 11554258485616, 6499270398159, 3348108992991, 1575580702584, 675248872536, 262596783764, 92263734836, 29135916264, 8217822536, 2054455634, 450978066, 85900584, 13983816, 1906884, 211876, 18424, 1176, 49, 1, 1};
            double ratiomax=0.,ratio;

            if(it_nbhood == 3) {
                binomialp = binomial9;
            }
            else if (it_nbhood == 5) {
                binomialp=binomial25;
            }
            else {
                binomialp=binomial49;
            }

            get_gliderinfo(gliderinfo, it_nbhood2*8);

            for (int i=0; i<it_nbhood2*8; i++) {
                if((i%it_nbhood2) == it_nbhood2-1) {
                    gliderinfo[i]=0;
                    ratio = 0.;
                }
                else ratio = (double) gliderinfo[i] / (double) binomialp[i%it_nbhood2];
                if(ratio > ratiomax) {
                    ratiomax = ratio;
                }
            }
            
            // fprintf(stderr,"ratiomax %g\n",ratiomax);
            
            for (ij=N2>>1;ij<N2;ij++) cgolg[ij] = 0;
            for (int i=0; i<it_nbhood2*8; i++) {
                ratio = (double) gliderinfo[i] / (double) binomialp[i%it_nbhood2];
                if((i%it_nbhood2) == it_nbhood2-1) {
                    for (int j=0;j< (N>>1);j++)
                        cgolg[((N-j)<<log2N)+i] = 0x7f7f7fff;
                }
                else {
                    for (int j=0;  j< (int) ((ratio*(double)(N>>1))/ratiomax)  ;j++)
                        cgolg[((N-j)<<log2N)+i] = 0xff00ff + (((i*0x1f/it_nbhood2))<<24)-(((i*0x1f/it_nbhood2))<<16);
                }
            }
            
        }
    }
    else if(colorfunction==9) {                                     // colorfunction based on unique labelling of separate components in image
      if (diagnostics & diag_component_labels) {
        for (ij=0; ij<NN2; ij++) labelcc[ij]=relabel[label[ij]];    // transfer labels to interactive working area
        if(xdisplay>=0 && ydisplay>=0) {
            if ((labelxy=label[xdisplay+ydisplay*N])) {
                for (ij=0; ij<NN2; ij++) if (label[ij]==complist[labelxy].label) labelcc[ij]=0xffff;    // label chosen component white at its current location
                d=complist[labelxy].log2n;
                short unsigned int conn = connlists[labelxy];
                while(conn) {
                    for (ij=0; ij<NN2; ij++) if (oldlabel[ij]==connections[conn].oldlab) {
                        if (labelcc[ij]==0xffff) labelcc[ij]=0xfffd; // color old connected components overlapping pink
                        else labelcc[ij]=0xfffe;                     // color old connected components not overlapping red
                    }
                    conn=connections[conn].next;
                }
                for (ij=0; ij<(1<<(d<<1)); ij++) labelcc[(ij&((1<<d)-1)) + (ij>>d)*N] = 0;  // initialize component drawing area to zero in top left corner
                if (complist[labelxy].quad) labelimage(complist[labelxy].quad, labelcc, 0xffff, 0); // extract this component and label it white
                else labelimage(complist[labelxy].patt, labelcc, 0xffff, 0);  // component involves patt not quad
            }
        }
 
        if (!gcolors) {
          for (ij=0; ij<NN2; ij++) {
            if (labelcc[ij]) {
                // quad = (uint64_t) relabel[label[ij]];            // use gene variable simply as label for connected component (no connection with gene)
                quad = (uint64_t) labelcc[ij];                      // use gene variable simply as label for connected component (no connection with gene)
                mask = quad * 11400714819323198549ull;              // map label (quasi-uniformly) into larger unsigned colour space
                mask = mask >> (64 - 32);                           // hash with optimal prime multiplicator down to 32 bits
                mask |= 0x080808ffull;                              // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                if (noveltyfilter && (diagnostics & diag_hash_patterns)) {
                    quad = complist[label[ij]].quad;
                    if((d=complist[label[ij]].patt)) popcount=smallpatts[d].activity; // number of times small (<= 4x4) pattern encountered previously
                    else if((q = (quadnode *) hashtable_find(&quadtable, quad)) != NULL) // if we reach here, quad should have been stored in hash table
                        popcount=q->activity;                       // number of times large pattern encountered previously (poss. also as part of larger patt)
                    else popcount=1;                                // should never occur, but just in case, assume novel
                    if (popcount>1) mask &= (mask&0x3f3f3fff);      // darken non novel components (currently a little too much)
                }
                if (labelcc[ij] == 0xffff) mask = 0xffffffffull;      // recolor selected component white
                else if (labelcc[ij] == 0xfffd) mask = 0xc0c0ffffull; // recolor connected old components overlapping with selected component pink
                else if (labelcc[ij] == 0xfffe) mask = 0x0000ffffull; // recolor connected old components not overlapping with selected component red
                cgolg[ij] = (int) mask;
            }
            else cgolg[ij] = 0;
          }
        }
        else {                                                       // with gcolors set, color according to computed inherited label colours
          for (ij=0; ij<NN2; ij++) {
            if (label[ij]) cgolg[ij] = (int) rgba(complist[label[ij]].gcolor);
            else cgolg[ij] = 0;
          }
        }
      } // diag_component_labels
    }
    else if(colorfunction==10){                                     //activities for patterns with size weighted colours
        uint64_t qid;
        for (ij=0; ij<NN2; ij++) labelcc[ij]=0;
        if(xdisplay>=0 && ydisplay>=0) {
            if ((qid=acttraceq[xdisplay+ydisplay*N])) {
                labelimage(qid, labelcc, 0xffff, 0);                // extract this component and label it white
            }
        }
        for (ij=0; ij<N2; ij++) {
            quad=acttraceq[ij];
            if(labelcc[ij]) mask = 0xffffffff;
            else if (quad == rootgene) mask = 0x3f3f3fff;           // grey color for background, all root genes
            else {
                if (activity_size_colormode == 0) {
                    mask = quad * 11400714819323198549ul;
                    mask = mask >> (64 - 32);                       // hash with optimal prime multiplicator down to 32 bits
                    mask |= 0x080808ffull;                          // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                }
                else if (activity_size_colormode == 1) {            // color by log2n, enclosing square size
                    if(acttraceqt[ij] && (q = (quadnode *) hashtable_find(&quadtable, quad)) != NULL) d = log2upper((unsigned int) q->size);
                    else if (!acttraceqt[ij] && quad<65536ull && smallpatts[quad].activity) {
                        d = log2size((short unsigned int) quad);
                    }
                    else {fprintf(stderr,"quad pattern not found in colorfunction for activities\n");d=0;}
                    if(d>log2N){fprintf(stderr,"Error in colorfunction 10, size error %d\n",d);d=0;}
                    setcolor(color,d);
                    mask = (color[0]<<8) | (color[1]<<16) |  (color[2]<<24) | 0xff;
                }
                else {                                              // color by sqrt of nr of live pixels (up to max value of 255)
                    int popmax = 255;
                    if(acttraceqt[ij] && (q = (quadnode *) hashtable_find(&quadtable, quad)) != NULL) popcount = q->pop1s;
                    else if (!acttraceqt[ij] && quad<65536ull && smallpatts[quad].activity) {POPCOUNT64C((uint64_t) quad,popcount);}
                    else {fprintf(stderr,"quad pattern not found in colorfunction for activities\n");popcount=0;}
                    if (activity_size_colormode == 3) popcount = (int) sqrtupper(popcount);
                    // fprintf(stderr,"step %d ij %d popcount %d\n",totsteps,ij,popcount);
                    if(popcount>popmax) popcount=popmax;
                    color[0]=popcount;color[0]= color[0]>255 ? 255: color[0];
                    color[1]=popcount<<3;color[1]=color[1]>255 ? 255: color[1];
                    color[2]=popcount<<6;color[2]=color[2]>255 ? 255: color[2];
                    for(d=0,mask=0xff;d<3;d++) mask |= color[d]<<((d<<3)+8);
                }
            }
            if ((npopulation[ij&Nmask]>>log2N)==(ij>>log2N)) mask = 0xffffffff;  // overlay plot with trace of density in white (except if pop=N2)
            cgolg[ij]= (int) mask;
        }
    }
    else if(colorfunction==11){                                     //genealogy based colours of ancestors at genealogycoldepth
    
        uint64_t ancestor,root;
        int ancestortypec;
        if((colorfunction2==-1) || (ancestortype !=2)) ancestortypec = ancestortype;
        else if (winnr) ancestortypec = 1;                          // ancestortype set to 3: use values 0 and 2 in two windows
        else ancestortypec = 0;
        if (ancestortypec==1) root = rootclone;
        else root = rootgene;

        for (ij=0; ij<N2; ij++) {
            if (gol[ij] && (diagnostics & diag_hash_genes)) {
                gene = (ancestortypec==2) ? ggolb[ij] : ggolg[ij];    // variable gene holds either birthid (clones) or gene at ij;
                ancestor=gene;
                for (j=1;j<=genealogycoldepth;j++) {
                    if((ancestortypec!=1) && (ancestor==rootgene)) break;                       // reached root, exit j loop
                    if((ancestortypec==1) && (ancestor&rootclone)) break;                       // reached root, exit j loop
                    else {
                        gene = ancestor;
                        if (ancestortypec==1) {
                            if((clonedataptr = (clonedata *) hashtable_find(&clonetable, gene)) != NULL) {
                                ancestor=clonedataptr->parentid;
                            }
                            else fprintf(stderr,"ancestor not found in clonealogies\n");
                        }
                        else {
                            if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
                                ancestor=genedataptr->firstancestor;
                            }
                            else fprintf(stderr,"ancestor not found in genealogies\n");
                        }
                    }
                }
                if (gene == 0ull) gene = 11778ull;                  // random color for gene==0
                mask = gene * 11400714819323198549ull;
                mask = mask >> (64 - 32);                           // hash with optimal prime multiplicator down to 32 bits
                mask |= 0x080808ffull;                              // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                if((ancestortypec!=1) && (ancestor==rootgene)) mask = 0x3f3f3fff;                // grey color for rootclone
                if((ancestortypec==1) && (ancestor&rootclone)) mask = 0x3f3f3fff;                // grey color for rootclone

                cgolg[ij] = (int) mask;
            }
            else cgolg[ij] = 0;
        }
    }
    else if(colorfunction==12){                                     // colouring based on periodicity of dynamic record of 16 last states in golr for live genes
        short unsigned int dscale[16] = {0xff,0xcf,0x7f,0x4f,0x2f,0x27,0x1f,0x1d,0x1b,0x19,0x17,0x15,0x14,0x13,0x12,0x11};
        int dx,dy;
        if (firstdebug) {
            golr_digest (0x1111111111111111ull, &d0, &d1, &jper, &dx, &dy);fprintf(stderr,"golr digest %llx %d %d %d %d %d\n",0x1111111111111111ull, d0, d1, jper, dx, dy);
            golr_digest (0x1111111111111119ull, &d0, &d1, &jper, &dx, &dy);fprintf(stderr,"golr digest %llx %d %d %d %d %d\n",0x1111111111111119ull, d0, d1, jper, dx, dy);
            golr_digest (0xaec8aec8aec8aec8ull, &d0, &d1, &jper, &dx, &dy);fprintf(stderr,"golr digest %llx %d %d %d %d %d\n",0xaec8aec8aec8aec8ull, d0, d1, jper, dx, dy);
            firstdebug=0;
        }
        for (ij=0; ij<N2; ij++) {
            if (ggol[ij] && (diagnostics & diag_hash_genes)) {
                gdiff = gene = ggolr[ij];                            // variable gene holds dynamical record stored in golr
                golr_digest (gene, &d0, &d1, &jper, &dx, &dy);
                if (d0 > 7) mask = 0x080808ffull;                   // dark grey color for no significant periodic match found
                else {
                    mask = 0xffull;
                    if ((d0==d1) && (~gene&0x8ull)) {               // survival with constant mismatch: ie static pattern : colour dark red
                        mask |= 0x3f<<8;
                    }
                    else if (abs(dx)< 0.25*(16-jper) && abs(dy) < 0.25*(16-jper)) {
                        mask |= 0x7f<<8;                            // stationary pattern: slightly less dark red
                    }
                    else {
                        mask |= (jper*(8-d0))<<25;
                        mask |= ((15-jper)*(8-d0))<<17;
                        mask |= dscale[(d0<<1)+1]<<8;
                    }
                }
                cgolg[ij] = (int) mask;
            }
            else cgolg[ij] = 0;
        }
    }
    for (ij=0; ij<N2; ij++) {                                       // convert BGRA format (pygame) to ARGB (PySDL2)
        uint32_t c = cgolg[ij];
        cgolg[ij]=((c&0xff)<<24) | ((c&0xff00)<<8) | ((c&0xff0000)>>8) | ((c&0xff000000)>>24);
    }
}
//------------------------------------------------------------- selectone -------------------------------------------------------------------------------
extern inline void selectone_of_2(int s, uint64_t nb2i, int nb[], uint64_t golg[], uint64_t golb[],uint64_t * birth, uint64_t *newgene, uint64_t *parentid, unsigned int *kch) {
// birth is returned 1 if ancestors satisfy selection condition. Selection of which of two genes to copy is newgene. Non-random result.
    unsigned int k,d0,d1,d2,d3,dd,swap;                  // number of ones in various gene combinations
    uint64_t livegenes[2],gdiff,gdiff0,gdiff1;           // various gene combinations
    uint64_t gene2centre;                                // gene function centres in sequence space
    int g0011,g0110,prey,prey2,ijanc[2];

    for(k=0;k<2;k++) livegenes[k] = golg[ijanc[k]=nb[(nb2i>>(k<<2))&0x7]];
    POPCOUNT64C(livegenes[0],d0);
    POPCOUNT64C(livegenes[1],d1);
    switch (selection) {
        case 0:                                          // integer value of sequence as fitness
            *birth = (livegenes[0]^livegenes[1]) ? 1ull: 0ull; // birth condition is two genes different
            *newgene = livegenes[0]>livegenes[1] ?  livegenes[0] : livegenes[1]; // choose one with larger gene to replicate
            *parentid = (livegenes[0]>livegenes[1]) ?  golb[ijanc[0]] : golb[ijanc[1]];
            break;

        case 1:                                          // number of ones in sequence as fitness
            *birth = (d0^d1) ? 1ull: 0ull;               // birth condition is two genes different in number of ones
            *newgene = (d0>d1) ? livegenes[0] : livegenes[1];
            *parentid = (d0>d1) ? golb[ijanc[0]] : golb[ijanc[1]];
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
                *parentid = (d0==0 && d1==3) ? golb[ijanc[swap^0]] : golb[ijanc[swap^1]];
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
                *parentid = (d0==0 && d1==3) ? golb[ijanc[swap^0]] : golb[ijanc[swap^1]];
            }
            else *newgene = 0ull;
            break;
        case 4:                                          // birth if 2 genes cooperate : closer to all 0, all 1 targets than ncoding and closer to each other than 64-ncoding
            gdiff=livegenes[0]^livegenes[1];
            POPCOUNT64C(gdiff,dd);
            *birth = dd>0 && dd<(64-ncoding) && ((d0<ncoding && d1>64-ncoding) || (d1<ncoding && d0>64-ncoding))  ? 1ull: 0ull; // birth if 2 genes close enough to targets
            if (d0<ncoding) {if(d0>64-d1) swap=1;else swap=0;}                  // need 64-d1 != d0 to avoid asymmetry in direction               && (64-d1 != d0)
            else {if(64-d0>=d1) swap=1; else swap=0;}
            *newgene = livegenes[swap];
            *parentid = golb[ijanc[swap]];
            break;
        case 5:                                          // predator prey model: prey-prey evolves to all 0, predator to complement of prey
            gdiff=livegenes[0]^livegenes[1];
            gdiff1=livegenes[0]^(~livegenes[1]);
            POPCOUNT64C(gdiff1,dd);
            prey = (d0<32) || (d1<32);                   // prey present : newgene is one with less ones, 1 prey : predator wins
            prey2 = (d0<32) && (d1<32);                  // 2 prey : newgene is one with less ones, 1 prey : predator wins
            // *birth = (gdiff && prey && dd<ncoding) ? 1ull: 0ull;  // birth if different and >=1 prey and close enough match)
            *birth = (gdiff && prey2) || (prey && (!prey2) && (dd<ncoding)) ? 1ull: 0ull; // birth if different and >=1 prey and close enough match)
            *newgene = (prey2 ? ((d0<d1) ? livegenes[0] : livegenes[1]) : ((d0<32) ? livegenes[1] : livegenes[0]));
            *parentid = (prey2 ? ((d0<d1) ? golb[ijanc[0]] : golb[ijanc[1]]) : ((d0<32) ? golb[ijanc[1]] : golb[ijanc[0]]));
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
            *newgene = (g0011 && (d0!=d3)) ? ((d0<d3) ? livegenes[0] : livegenes[1]) : ((d2<d1) ? livegenes[0] : livegenes[1]);
            *parentid = (g0011 && (d0!=d3)) ? ((d0<d3) ? golb[ijanc[0]] : golb[ijanc[1]]) : ((d2<d1) ? golb[ijanc[0]] : golb[ijanc[1]]);
            break;
        case 7:                                          // neutral selection but selective birth only occurs if two chosen sequences are same (different) (NB uses RNG)
                                                         // note that birth is already suppressed for 3 identical live nbs unless enforcebirth-3 bit on
                                                         // we want to use this routine only for two live nbs and only when genes not same
            *birth = livegenes[0]^livegenes[1] ? 1ull: 0ull;
            *newgene = golg[nb[*kch]];
            *parentid = golb[nb[*kch]];
            break;
        // cases 8+: this subroutine is not used
        default:
            fprintf(stderr,"Error: two live gene fitness value %d is not implemented\n",selection);
            exit(1);
    }
}
//.......................................................................................................................................................
extern inline int selectone_of_s(unsigned int *kch, int s, uint64_t nb1i, int nb[], uint64_t golg[], uint64_t golb[], uint64_t golr[], uint64_t *birth, uint64_t *newgene, uint64_t *parentid, uint64_t *nbmask, int ij) {
// result is number of equally fit best neighbours that could be the ancestor (0 if no birth allowed, 1 if unique) and list of these neighbour indices
// birth is returned 1 if ancestors satisfy selection condition. Selection of which of genes to copy is newgene. Non-random result.
    unsigned int k,nbest,ijanc[8],kchs[8];                // index for neighbours and number in best fitness category and ij, kch for up to 8 possible ancestors
    unsigned int d[8],dmax[8],p[8],dS,dB,d0,d1,d2;         // d number of ones or mismatches, dmx max mismatches,p period, + other variables
    int psx[8],psy[8];                                    // periodic shift in x and y for optimal period of recorded displacements
    unsigned int scores[8];                               // cumulative scores for pairwise games of individual livegenes (used case repselect == 7)
    uint64_t livegenes[8],gdiff,extremval,bestnbmask,birthid;
    unsigned int repselect = (repscheme & R_47_repselect)>>4; //

    birthid = (uint64_t) totsteps;
    birthid = (birthid << 32)+rootclone+ij;
    if(s==0) {*birth = 1ull;  *newgene = genegol[selection-8]; *parentid = birthid; *kch = 0; return(1);}

    for(k=0;k<s;k++) {
        kchs[k]=(nb1i>>(k<<2))&0x7;
        ijanc[k] = nb[kchs[k]];
        livegenes[k] = golg[ijanc[k]];
        if(repselect<8) {POPCOUNT64C(livegenes[k],d[k]);}
        else golr_digest (golr[ijanc[k]], d+k, dmax+k, p+k, psx+k, psy+k);
    }

    switch (repselect) {
        case 0:                                          // integer value of sequence as fitness : smallest (case 0) or largest (case 1)
        case 1:
            if(repselect&0x1) for(extremval= 0ull,k=0;k<s;k++) extremval = (livegenes[k] >= extremval ? livegenes[k] : extremval); // find value of fittest gene max value
            else              for(extremval=~0ull,k=0;k<s;k++) extremval = (livegenes[k] <= extremval ? livegenes[k] : extremval); // find value of fittest gene min value
            for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= ((livegenes[k]==extremval) ? 1ull<< (k+0*nbest++) : 0ull); // find set of genes with equal best value
            *birth = (nbest>0) ? 1ull: 0ull;             // birth condition may include later that genes not all same
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;   // execute loop until first optimal live neighbour found at k (a one)
            if (k==s) {k=0;*birth = 0ull;}               // in case no genes with best value, no birth, avoid k being out of bounds below
            *newgene = livegenes[k&0x7];                 // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            break;
        case 2:
        case 3:                                          // number of ones in sequence as fitness : smallest (case 2) or largest (case 3)
            if(repselect&0x1) for(extremval= 0ull,k=0;k<s;k++) extremval = (d[k]>= extremval ? d[k] : extremval); // find value of fittest gene max nr ones
            else              for(extremval=~0ull,k=0;k<s;k++) extremval = (d[k]<= extremval ? d[k] : extremval); // find value of fittest gene min nr ones
            for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= (d[k]==extremval ? 1ull<<(k+0*nbest++) : 0ull); // find set of genes with equal best value
            *birth = ((nbest>0) ? 1ull: 0ull);           // birth condition may include later that genes not all same
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            if (k==s) {k=0;*birth = 0ull;}               // in case no genes with best value, no birth, avoid k being out of bounds below
            *newgene = livegenes[k&0x7];                 // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            break;
        case 4:                                          // neutral selection
            for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= 1ull<<(k+0*nbest++); // find set of genes with equal best value : in this case all of them
            *birth = 1ull;                               // birth condition is always true
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            if (k==s) {k=0;*birth = 0ull;}               // in case no genes with best value, no birth, avoid k being out of bounds below
            *newgene = livegenes[k&0x7];                 // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            break;
        case 5:                                          // neutral selection but birth only occurs if some sequences are different
            for(gdiff=0ull,k=1;k<s;k++) if ((gdiff=livegenes[0]^livegenes[k])) break; // test whether all genes the same
            if (gdiff) for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= 1ull<<(k+0*nbest++); // find set of genes with equal best value : in this case all of them
            else { nbest = 0; bestnbmask = 0ull;}
            *birth = gdiff ? 1ull : 0ull;                // birth condition is genes not all same
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            if (k==s) {k=0;*birth = 0ull;}               // in case no genes with best value, no birth, avoid k being out of bounds below
            *newgene = livegenes[k&0x7];                 // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            break;
        case 6:                                          // penalize genes by a cost for the numebr of LUT entries realized for survival and/or birth
          switch(selection) {
            case 8:                                      // totalistic lut penalty of gene in fixed length encoding : first 8 bits survival, next 8 bits birth
                for(k=0;k<s;k++) {
                    // livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    POPCOUNT64C(livegenes[k]&0xffull,dS);
                    POPCOUNT64C((livegenes[k]>>8)&0xffull,dB);
                    d[k]=24-dS-(dB<<1);
                }
                break;
            case 9:                                      // totalistic lut penalty of gene in variable length encoding : 4 bit patterns S 0xxx  B 1xxx
                for(k=0;k<s;k++) {
                    //livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    PATTERN4(livegenes[k],0x0,d0)
                    PATTERN4(~livegenes[k]&0x8888888888888888ull,0x8,d1)
                    PATTERN4(livegenes[k]&0x8888888888888888ull,0x8,d2)
                    d[k]=48-(d1-d0)-2*d2;
                }
                break;
            case 10:
                for(k=0;k<s;k++) {
                    // livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    POPCOUNT64C(livegenes[k]&0x7fffffull,dS);        // 23 coding bits for survival
                    POPCOUNT64C((livegenes[k]>>32)&0x7fffffull,dB);  // 23 coding bits for birth
                    d[k]=57-dS-(dB<<1);
                }
                break;
            case 11:                                     // dist. lut penalty of gene in variable length encoding : 8 bit patterns S 0xxxrrrr  B 1xxxrrrr
                for(k=0;k<s;k++) {
                    // livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    PATTERN8(livegenes[k],0x00,d0)
                    PATTERN8(~livegenes[k]&0x8080808080808080ull,0x80,d1)
                    PATTERN8(livegenes[k]&0x8080808080808080ull,0x80,d2)
                    d[k]=24-(d1-d0)-2*d2;
                }
                break;
            case 12:                                     // canon. lut penalty in fixed length encoding 32 bit survival 32 bit birth
            case 14:
                 for(k=0;k<s;k++) {
                    livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    POPCOUNT64C(livegenes[k]&0xffffffffull,dS);        // 32 coding bits for survival
                    POPCOUNT64C((livegenes[k]>>32)&0xffffffffull,dB);  // 32 coding bits for birth
                    d[k]=96-dS-(dB<<1);
                }
                break;
            case 13:                                    // canon. lut penalty of gene in variable length encoding : 8 bit patterns S 0xxxwrrr  B 1xxxwrrr
                for(k=0;k<s;k++) {                      // could be made more accurate by accounting for upper pair bit states for rest of 5-subset too
                    // livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    PATTERN8(livegenes[k]&0xffffffffffff,0x00,d0)
                    PATTERN8(~livegenes[k]&0x808080808080ull,0x80,d1)
                    PATTERN8(livegenes[k]&0x808080808080ull,0x80,d2)
                    d[k]=24-(d1-d0)-2*d2;
                }
                break;
            case 15:                                    // dist. or canon. lut penalty of gene in variable length encoding : 8 bit patterns S 0xxxwwrr  B 1xxxwwrr
                for(k=0;k<s;k++) {                      // could be made more accurate by accounting for upper quartet bit states for rest of 6-subset too
                    // livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    PATTERN8(livegenes[k]&0xffffffffff,0x00,d0)
                    PATTERN8(~livegenes[k]&0x8080808080ull,0x80,d1)
                    PATTERN8(livegenes[k]&0x8080808080ull,0x80,d2)
                    d[k]=24-(d1-d0)-2*d2;
                }
                break;
            default:
                fprintf(stderr,"Error: selectone of s for selection %d is not implemented\n",selection);
                exit(1);
          }
          for(extremval=0ull,k=0;k<s;k++) extremval = d[k]>= extremval ? d[k] : extremval; // find value of fittest gene
          for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= (d[k]==extremval) ? 1ull<<(k+0*nbest++) : 0ull; // find set of genes with equal best value
          *birth = 1ull;                                // birth condition is always met since extremval set is always > 0
          for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
          if (k==s) {k=0;*birth = 0ull;}                // in case no genes with best value, no birth, avoid k being out of bounds below
          // if(k==s) fprintf(stderr,"Error in selectone of s case 6: k>=s (%d > %d)\n",k,s);
          *newgene = livegenes[k&0x7];                  // choose first of selected set to replicate (can make positional dependent choice instead externally)
          *parentid=golb[ijanc[k&0x7]];
          *kch = kchs[k];
          break;
        case 7:                                         // scissors-stone-well-paper game on number ones mod 4, scissors 0 stone 1 well 2 paper 3
            for (k=0;k<s;k++) {                         // exception to numerical order: sc>pa. Do tournament with all others for each livegene and score
                scores[k]=0; d0 = d[k]&0x3;
                for (int k1=0;k1<s;k1++) {
                    d1 = d[k1]&0x3;                     // if k==k1 there is zero score in next line so allow it
                    scores[k] += ( (d0 > d1 && !(d0 == 3 && d1 == 0)) || (d0 == 0 && d1 == 3) ) ? 1 : 0;
                }
            }
            for(extremval= 0ull,k=0;k<s;k++) extremval = (scores[k]>= extremval ? scores[k] : extremval); // find value of fittest gene, best score in tournament
            if (extremval) {
                for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= (scores[k]==extremval ? 1ull<<(k+0*nbest++) : 0ull); // find set of genes with equal best value
            }
            else {   // no positive scores, no birth
                bestnbmask=0ull;nbest=0;
            }
            *birth = nbest ? 1ull: 0ull;               // birth condition
            if (s==3) {                                // included for compatibility with update_23() for enforcebirth off
                for(d0=k=0;k<s;k++) for (int k1=0;k1<k;k1++) d0 |= livegenes[k]^livegenes[k1] ? 1 : 0; //  d0 = 1 if  genes are not all the same
                *birth = d0 ? 1ull: 0ull;              // birth condition modified to only be true if some genes difft
                nbest = d0 ? 3 : 0;
                bestnbmask = d0 ? 0x7ull : 0;
            }
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            if (k==s) {k=0;*birth = 0ull;}             // in case no genes with best value, no birth, avoid k being out of bounds below

            *newgene = livegenes[k&0x7];               // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            //if(ij==IJDEBUG) fprintf(stderr,"DEBUG In selectone_of_s at totsteps=%d ij=%d s=%d nbest=%d extremval=%llx bestnbmask=%llx nb1i=%llx kch=%d\n",totsteps,ij,s,nbest,extremval,bestnbmask,nb1i,*kch);
            break;
        case 8:                                        // case 8-15 are intended for golr selection modes
        case 9:
            if(repselect&0x1) for(extremval= 0ull,k=0;k<s;k++) extremval = (p[k]>= extremval ? p[k] : extremval); // find value of fittest gene max period
            else              for(extremval=~0ull,k=0;k<s;k++) extremval = (p[k]<= extremval ? p[k] : extremval); // find value of fittest gene min period
            for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= (p[k]==extremval ? 1ull<<(k+0*nbest++) : 0ull); // find set of genes with equal best value
            *birth = ((nbest>0) ? 1ull: 0ull);           // birth condition may include later that genes not all same
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            if (k==s) {k=0;*birth = 0ull;}               // in case no genes with best value, no birth, avoid k being out of bounds below
            *newgene = livegenes[k&0x7];                 // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            break;
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
        case 15:
            if(repselect&0x1) for(extremval= 0ull,k=0;k<s;k++) extremval = (d[k]>= extremval ? d[k] : extremval); // find value of fittest gene most aperiodic
            else              for(extremval=~0ull,k=0;k<s;k++) extremval = (d[k]<= extremval ? d[k] : extremval); // find value of fittest gene precise period
            for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= (d[k]==extremval ? 1ull<<(k+0*nbest++) : 0ull); // find set of genes with equal best value
            *birth = ((nbest>0) ? 1ull: 0ull);           // birth condition may include later that genes not all same
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            if (k==s) {k=0;*birth = 0ull;}               // in case no genes with best value, no birth, avoid k being out of bounds below
            *newgene = livegenes[k&0x7];                 // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            break;
        default:
            fprintf(stderr,"Error: s = %d live gene repselect %d is not implemented\n",s,repselect);
            exit(1);
    }
    for (*nbmask=0ull,k=0;k<s;k++)
        *nbmask |= ((bestnbmask>>k)&0x1ull)<<((nb1i>>(k<<2))&0x7);

    return(nbest);
}
//------------------------------------------------------------- selectone_nbs ---------------------------------------------------------------------------
extern inline void selectone_nbs(int s, uint64_t nb2i, int nb[], uint64_t gol[], uint64_t golg[],uint64_t golb[], uint64_t * birth, uint64_t *newgene, uint64_t *parentid, unsigned int *kch) {
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
        k = (nb2i>>(kanc<<2))&0x7;
        ij1 = nb[k];
        *parentid=golb[ij1];
        *newgene = golg[ij1];
        *kch = k;
    }
    // else no change to birth or kch
}
//-------------------------------------------------------------- selectdifftx ---------------------------------------------------------------------------
extern inline unsigned int selectdifft0(uint64_t nbmask, int *crot, int *kodd) {
    *kodd = 0;
    *crot = 0;
    return(0);                                                 // replication of live nb in bit 0 of canonical rotation
}
//.......................................................................................................................................................
extern inline unsigned int selectdifft1(uint64_t nbmask, int *crot, int *kodd) {
// selection based on canonical rotation
    int k,kmin;
    uint64_t nbmaskr, nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {         // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);      // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                              // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                // neighbor mask rotate min is current rotation
            kmin = k;                                          // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;
    *crot = 0;
    return(kmin);                                              // replication of live nb in bit 0 of canonical rotation
}
//.......................................................................................................................................................
extern inline unsigned int selectdifft2(uint64_t nbmask, int *crot, int *kodd) {
// selection based on canonical rotation to bunched pair, choose clockwise or anti-clockwise one (for R_2_canonical_nb)
    int k,kmin;
    uint64_t nbmaskr, nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {         // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);      // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                              // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                // neighbor mask rotate min is current rotation
            kmin = k;                                          // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;                                        // if not canonical, replication of live neighbour in other (non zero k) position
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

    if (canonical) return(kmin);                               // replication of live neigbour in bit 0 of canonical rotation
    else return((kmin+k)&0x7);                                 // rotate unique nb k left (kmin) back to orig nb pat
}
//.......................................................................................................................................................
extern inline unsigned int selectdifft3(uint64_t nbmask, int *crot, int *kodd) {
    unsigned int k,kmin;
    uint64_t nbmaskr,nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {         // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 1ull)<<7) | (nbmaskr>>1);        // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                              // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                // neighbor mask rotate min is current rotation
            kmin = k;                                          // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;                                        // replication of live neighbour in most different position
    switch (nbmaskrm) {                                        //              x07    x0b    x0d    x13    x15    x19    x25
        case 0x07ull : k = 1; *crot = 0; break;                // 00000111    |012|  <-
        case 0x0bull : k = 0; *crot = 1; break;                // 00001011    |...|  |01.|  <-
        case 0x0dull : k = 3; *crot = 2; break;                // 00001101    |...|  |..3|  |0.2|   <-
        case 0x13ull : k = 1; *crot = 3; break;                // 00010011           |...|  |..3|  |01.|   <-
        case 0x15ull : k = 2; *crot = 4; break;                // 00010101                  |...|  |...|  |0.2|   <-
        case 0x19ull : k = 0; *crot = 5; break;                // 00011001                         |..4|  |...|  |0..|   <-
        case 0x25ull : k = 5; *crot = 6; break;                // 00100101                                |..4|  |..3|  |0.2|  <-
        default  : {                                           //                                                |..4|  |...|
                                                               //                                                       |.5.|
            fprintf(stderr,"Error in canonical rotation for three live neighbours nbmaskrm = %llx for mask %llx\n",nbmaskrm,nbmask); k = 0;
        } //default case
    } //switch
    if (canonical) return(kmin);                               // replication of live neigbour in bit 0 of canonical rotation
    else return((kmin+k)&0x7);                                 // rotate unique nb k left (kmin) back to orig nb pat

}
//.......................................................................................................................................................
extern inline unsigned int selectdifft4(uint64_t nbmask, int *crot, int *kodd) {
    int k,kmin;
    uint64_t nbmaskr,nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {         // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);      // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                              // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                // neighbor mask rotate min is current rotation
            kmin = k;                                          // no of times rotated to right
        }
    }
    *kodd = kmin&0x1;                                          // replication of live neighbour in most central position (left disambiguation)
    switch (nbmaskrm) {                                        //              x07f    x17    x1b    x1d    x27    x2b    x2d   x33   x35    x55
        case 0x0full : k = 1; *crot = 0; break;                // 00001111    |012|  <-
        case 0x17ull : k = 2; *crot = 1; break;                // 00010111    |..3|  |012|  <-
        case 0x1bull : k = 1; *crot = 2; break;                // 00011011    |...|  |...|  |01.|   <-
        case 0x1dull : k = 2; *crot = 3; break;                // 00011101           |..4|  |..3|  |0.2|   <-
        case 0x27ull : k = 2; *crot = 4; break;                // 00100111                  |..4|  |..3|  |012|   <-
        case 0x2bull : k = 3; *crot = 5; break;                // 00101011                         |..4|  |...|  |01.|   <-
        case 0x2dull : k = 2; *crot = 6; break;                // 00101101                                |.5.|  |..3|  |0.2|  <-
        case 0x33ull : k = 1; *crot = 7; break;                // 00110011                                       |.5.|  |..3|  |01.|  <-
        case 0x35ull : k = 2; *crot = 8; break;                // 00110101                                              |.5.|  |...|  |0.2|  <-
        case 0x55ull : k = 2; *crot = 9; break;                // 01010101                                                     |.54|  |...|  |0.2|  <-
        default  : {                                           //                                                                     |.54|  |...|
                                                               //                                                                            |6.4|
            fprintf(stderr,"Error in canonical rotation for four live neighbours nbmaskrm = %llx for mask %llx\n",nbmaskrm,nbmask); k = 0;
        } //default case
    } //switch
    if (canonical) return(kmin);                               // replication of live neigbour in bit 0 of canonical rotation
    else return((kmin+k)&0x7);                                 // rotate unique nb k left (kmin) back to orig nb pat
}
//.......................................................................................................................................................
extern inline unsigned int selectdifft5(uint64_t nbmask, int *crot, int *kodd) {
    unsigned int k,kmin;
    uint64_t nbmaskr,nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {         // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 1ull)<<7) | (nbmaskr>>1);        // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                              // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                // neighbor mask rotate min is current rotation
            kmin = k;                                          // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;
    switch (nbmaskrm) {                                        //              x1f    x2f    x3d    x3b    x57    x37    x5b
        case 0x1full : k = 2; *crot = 0; break;                // 00011111    |012|  <-
        case 0x2full : k = 2; *crot = 1; break;                // 00101111    |..3|  |012|  <-
        case 0x3dull : k = 3; *crot = 2; break;                // 00111101    |..4|  |..3|  |0.2|   <-
        case 0x3bull : k = 3; *crot = 3; break;                // 00111011           |.5.|  |..3|  |01.|   <-
        case 0x57ull : k = 2; *crot = 4; break;                // 01010111                  |.54|  |..3|  |012|   <-
        case 0x37ull : k = 2; *crot = 5; break;                // 00110111                         |.54|  |...|  |012|   <-
        case 0x5bull : k = 3; *crot = 6; break;                // 01011011                                |6.4|  |...|  |01.|  <-
        default  : {                                           //                                                |.54|  |..3|
                                                               //                                                       |6.4|
            fprintf(stderr,"Error in canonical rotation for five live neighbours nbmaskrm = %llx for mask %llx\n",nbmaskrm,nbmask); k = 0;
        } //default case
    } //switch
    if (canonical) return(kmin);                               // replication of live neigbour in bit 0 of canonical rotation
    else return((kmin+k)&0x7);                                 // rotate unique nb k left (kmin) back to orig nb pat
}
//.......................................................................................................................................................
extern inline unsigned int selectdifft6(uint64_t nbmask, int *crot, int *kodd) {
// selection based on canonical rotation to bunched pair, choose clockwise or anti-clockwise one (for R_2_canonical_nb)
    int k,kmin;
    uint64_t nbmaskr, nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {         // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);      // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                              // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                // neighbor mask rotate min is current rotation
            kmin = k;                                          // no of times rotated to right
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
    if (canonical) return(kmin);                               // replication of live neigbour in bit 0 of canonical rotation
    else return((kmin+k)&0x7);                                 // rotate unique nb k left (kmin) back to orig nb pat
}
//.......................................................................................................................................................
extern inline unsigned int selectdifft7(uint64_t nbmask, int *crot, int *kodd) {
// selection based on canonical rotation
    int k,kmin;
    uint64_t nbmaskr, nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {         // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);      // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                              // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                // neighbor mask rotate min is current rotation
            kmin = k;                                          // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;
    *crot = 0;
    if (canonical) return(kmin);                               // replication of live nb in bit 0 of canonical rotation
    else  return((kmin+3)&0x7);                                // replication of live nb in bit 3 (middle) of canonical rotation
}
//.......................................................................................................................................................
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
                             // if (*crot==2) *nsame = 2;      // this case can be resolved with good rotation symmetry
                             if (*crot==7) *nsame = 2;
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
//.......................................................................................................................................................
extern inline void analyze_nbs(int ij, uint64_t gol[], int nnb[], uint64_t *nnb1i, int *ns) {
    int i,j,ip1,im1,jp1,jm1,s,k;
    uint64_t nb1i,gols;
    
    i = ij & Nmask;  j = ij >> log2N;                                       // row & column
    jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                       // toroidal (j+1)*N and (j-1)*N
    ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                             // toroidal i+1, i-1
    nnb[0]=jm1+im1; nnb[1]=jm1+i; nnb[2]=jm1+ip1; nnb[3]=j*N+ip1;           // new order of nbs
    nnb[4]=jp1+ip1; nnb[5]=jp1+i; nnb[6]=jp1+im1; nnb[7]=j*N+im1;
    for (s=0,nb1i=0ull,k=0;k<8;k++) {                                       // packs non-zero nb indices in first up to 8*4 bits
        gols=gol[nnb[k]];                                                   // whether neighbor is alive
        s += gols;                                                          // s is number of live nbs
        nb1i = (nb1i << (gols<<2)) + (gols*k);                              // nb1i is packed list of live neighbour indices
    }
    *nnb1i = nb1i;
    *ns = s;
}
//.......................................................................................................................................................
extern inline uint64_t disambiguate(unsigned int *kchx, uint64_t nb1i, int nb[],  uint64_t gol[], uint64_t golg[], uint64_t golb[], int nsame, uint64_t *birth, uint64_t *parentid, uint64_t *ancestor, int ij) {
    uint64_t gene,newgene,randnr;
    int k,ijanc,nnb[8],ns,discase,deathlikely;
    unsigned int kch,kch1;
    uint64_t nnb1i;
    
    kch = *kchx;
    discase=(repscheme>>8)&0x7;
    switch (discase) {
        case 0:  RAND128P(randnr);                                                      // random choice
                 kch += ((nsame-1)&(randnr>>32))<< (nsame ==4 ? 1 : 2);
                 kch &= 0x7;*kchx = kch;
                 ijanc = nb[kch];
                 *parentid = golb[ijanc];
                 *ancestor = golg[ijanc];
                 return( golg[ijanc]);
        case 1:  ijanc = nb[kch];                                                        // ignore asymmetry issue, continue regardless;
                 *parentid = golb[ijanc];
                 // if(*parentid == 0ull) fprintf(stderr,"error in disambiguate case 1: parentid set to golb[%d] which is 0 kch %d\n",ijanc,kch);
                 *ancestor = golg[ijanc];
                 return( golg[ijanc]);
        case 2:  *birth = 0ull; return(0ull);                                            // abandom birth attempt
        case 3:  *parentid = (((uint64_t) totsteps) <<32) + rootclone + ij;              // choose one GoL input gene, default ancestor for input genes
                 *ancestor = rootgene;
                 return(genegol[selection-8]);
        case 4:  for (newgene=golg[nb[kch]],kch1=kch,k=1;k<nsame;k++) {                   // choose minimum value gene
                     kch1+=k*(nsame==4 ? 2 : 4);kch1 &= 0x7;
                     ijanc=nb[kch1];
                     gene=golg[ijanc];
                     if (gene<newgene) {newgene = gene;kch=kch1;}
                 };
                 *kchx = kch;
                 ijanc=nb[kch];
                 *parentid = golb[ijanc];
                 *ancestor = newgene;
                 return(newgene);
        case 5:
        case 6:  for (newgene=golg[nb[kch]],kch1=kch,k=0;k<nsame;k++) {                        // choose first gene which is likely to die
                     kch1+=k*(nsame==4 ? 2 : 4);
                     kch1 &= 0x7;
                     ijanc=nb[kch1];
                     analyze_nbs(ijanc, gol, nnb, &nnb1i, &ns);
                     if(discase==5) deathlikely = ns<2 || ns>3;
                     else           deathlikely = ns<2 || ns>3 || ((overwritemask>>(ns-1))&0x1ull);
                     if(deathlikely) {                           // approx. substitute for death/overwrite calculation
                        newgene = golg[ijanc];
                        kch=kch1;
                        break;
                     }
                 };
                 *kchx = kch;
                 ijanc=nb[kch];
                 *parentid = golb[ijanc];
                 *ancestor = newgene;
                 return(newgene);
        case 7:  *parentid = (uint64_t) totsteps; *parentid = (*parentid <<32) + rootclone + ij;// default ancestor for input genes
                 RAND128P(randnr);                                                      // random choice
                 // *kchx = kch;    no change, retains kch unaltered as in case 1
                 *ancestor = rootgene;
                 return(randnr);                                                         // choose random gene
        default: fprintf(stderr,"Error in switch of ambiguous rotation resolution, should never reach here\n");
                 return(0ull);
    }
}
//------------------------------------------------------------ hash gene inline fns ---------------------------------------------------------------------
extern inline void hashaddgene(int ij,uint64_t gene,uint64_t ancestor,uint64_t *golb,uint64_t parentid,uint64_t mutation) {
    genedata gdata;
    uint64_t birthid;
    extern inline void hashaddclone(uint64_t birthid, uint64_t parentid, uint64_t gene);
    
    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
        genedataptr->lasttime = totsteps;
        if(mutation) {
            genedataptr->recenttime = (short unsigned int) totsteps;
        }
        genedataptr->popcount++;
    }
    else {
        gdata=ginitdata;
        gdata.gene = gene;
        gdata.firsttime = gdata.recenttime = gdata.lasttime = (short unsigned int) totsteps;
        gdata.firstancestor = ancestor;
        hashtable_insert(&genetable, gene,(genedata *) &gdata);
    }
    if(ancestor != rootgene) {
        if((genedataptr = (genedata *) hashtable_find(&genetable, ancestor)) == NULL) {
            fprintf(stderr,"error in hashaddgene step %d ij %d, the ancestor %llx of gene %llx to be stored is not stored (mutation %llx)\n",totsteps,ij,ancestor,gene,mutation);
        }
    }
    if(diagnostics & diag_hash_clones) {
        if (mutation || (parentid&rootclone)) {                   // when s=0 birth involves the creation of a new clone from root and rootclone bit set in parentid
            birthid = ((uint64_t) totsteps)<<32;
            birthid |= ij;
            hashaddclone(birthid,parentid,gene);
            *golb = birthid;
        }
        else {
            *golb = parentid;
            if((clonedataptr = (clonedata *) hashtable_find(&clonetable, parentid)) != NULL) clonedataptr->popln++;
            else {
                fprintf(stderr,"step %d error in hashclone update, %llx clone not saved\n",totsteps,parentid);
                fprintf(stderr,"ij %llx rootclone %d tstep %lld\n",parentid&N2mask,(parentid&rootclone)?1:0,parentid>>32);
            }
        }
    }
}
//.......................................................................................................................................................
extern inline void hashdeletegene(uint64_t gene,uint64_t birthid,const char errorformat[]) {
    extern inline void hashdeletefromclone(uint64_t birthid);
    
    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {genedataptr->popcount--;}
    else fprintf(stderr,errorformat,totsteps,gene);     // errorformat must contain %d and %llx format codes in this order
    
    if(diagnostics & diag_hash_clones) hashdeletefromclone(birthid);
}
//.......................................................................................................................................................
extern inline void hashgeneextinction(uint64_t gene,const char errorformat[]) {
    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
        if(genedataptr->popcount == 0) {
            genedataptr->lastextinctiontime = (short int) totsteps;
            genedataptr->nextinctions++;
        }
    }
    else fprintf(stderr,errorformat,3,totsteps,gene);
}
//.......................................................................................................................................................
extern inline void hashgeneactivity(uint64_t gene, const char errorformat[]) {
        if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) genedataptr->activity++;
        else fprintf(stderr,errorformat,4,totsteps,gene);
}
//------------------------------------------------------------ hash clone inline fns --------------------------------------------------------------------
extern inline void hashaddclone(uint64_t birthid, uint64_t parentid, uint64_t gene) {
    clonedata cdata;

    if((clonedataptr = (clonedata *) hashtable_find(&clonetable, birthid)) != NULL) {
        fprintf(stderr,"error in hashaddclone, %llx already present\n",birthid);
    }
    else {
        cdata=cinitdata;
        cdata.birthid=birthid;
        cdata.parentid=parentid;
        cdata.gene = gene;
        hashtable_insert(&clonetable, birthid,(clonedata *) &cdata);
    }
}
//.......................................................................................................................................................
extern inline void hashdeletefromclone(uint64_t birthid) {
    if((clonedataptr = (clonedata *) hashtable_find(&clonetable, birthid)) != NULL) {
        if(clonedataptr->popln) clonedataptr->popln--;
        else fprintf(stderr,"error in deleting individual from clone: popln already zero\n");
        // if (!clonedataptr->popln) hashtable_remove( &clonetable, birthid );
    }
    else fprintf(stderr,"error deleting from clone with birthid %llx at time %d : not found in hash table\n",birthid,totsteps);
}
//.......................................................................................................................................................
extern inline void hashcloneactivity(uint64_t birthid, const char errorformat[]) {
        if((clonedataptr = (clonedata *) hashtable_find(&clonetable, birthid)) != NULL) clonedataptr->activity++;
        else fprintf(stderr,errorformat,6,totsteps,birthid);
}
//------------------------------------------------------- hash quadtree inline fns ----------------------------------------------------------------------
extern inline uint16_t rotate16(uint16_t patt) {                                        // rotate bits in 4x4 pattern for 90 deg clockwise rotation
    uint16_t rotate4[16]={0,2,8,10,1,3,9,11,4,6,12,14,5,7,13,15};
    uint16_t nw,ne,sw,se;
    
    nw=patt&0xf;patt>>=4; nw=rotate4[nw];
    ne=patt&0xf;patt>>=4; ne=rotate4[ne];
    sw=patt&0xf;patt>>=4; sw=rotate4[sw];
    se=patt&0xf;          se=rotate4[se];
    return sw | (nw<<4) | (se<<8) | (ne<<12) ;
}
//.......................................................................................................................................................
extern inline uint64_t rotate64(uint64_t patt) {                                        // rotate bits in 8x8 pattern for 90 deg clockwise rotation
    uint64_t nw,ne,sw,se;
    
    nw=patt&0xff;patt>>=16; nw=(uint64_t) rotate16((uint16_t) nw);
    ne=patt&0xff;patt>>=16; ne=(uint64_t) rotate16((uint16_t) ne);
    sw=patt&0xff;patt>>=16; sw=(uint64_t) rotate16((uint16_t) sw);
    se=patt&0xff;           se=(uint64_t) rotate16((uint16_t) se);
    return sw | (nw<<16) | (se<<32) | (ne<<48) ;
}
//.......................................................................................................................................................
extern inline void rotate4x64(uint64_t *nw, uint64_t *ne, uint64_t *sw, uint64_t *se) { // rotate bits in 16x16 pattern for 90 deg clockwise rotation
    *nw=rotate64(*sw); *ne=rotate64(*sw); *sw=rotate64(*se); *se=rotate64(*ne);
}
//.......................................................................................................................................................
extern inline void rotatequad(uint64_t *nw, uint64_t *ne, uint64_t *sw, uint64_t *se) { // rotate bits in quad pattern for 90 deg clockwise rotation
// NYI 1. hash_node(*nw,*ne,*sw,*se) 2. lookup hash entry of hash & check id 3. if patt then call rotate4x64 else do 4 calls to rotatequad with subnodes
}
//.......................................................................................................................................................
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
extern inline quadnode * hash_patt8_store(const uint64_t h, const uint64_t patt) {
    int nr1;
    quadinit.hashkey = h;
    quadinit.isnode=0;
    quadinit.nw=patt;
    quadinit.ne=0;
    quadinit.sw=0;
    quadinit.se=0;
    quadinit.size= 8;
    quadinit.firsttime=totsteps;
    quadinit.lasttime=totsteps;
    POPCOUNT64C(patt,nr1);
    quadinit.pop1s =nr1;
    hashtable_insert(&quadtable, h,(quadnode *) &quadinit);
    return (quadnode *) hashtable_find(&quadtable, h);
}
//.......................................................................................................................................................
extern inline void hash_patt4_find(const short unsigned int patt) {

    if(smallpatts[patt].activity) {                                         // pattern found
        smallpatts[patt].activity++;
        smallpatts[patt].lasttime=totsteps;
    }
    else {                                                                  // store new pattern
        // smallpatts[patt].size=log2upper(patt);
        smallpatts[patt].activity++;
        smallpatts[patt].firsttime=totsteps;
        smallpatts[patt].lasttime=totsteps;
    }
}
//.......................................................................................................................................................
extern inline quadnode * hash_patt8_find(const uint64_t patt) {
        quadnode *q;
        uint64_t h,npatt;
        short unsigned int nw,ne,sw,se;
        const uint64_t randomizer = 11400714819323198549ull;

        h = patt_hash(patt,patt+1,patt+2,patt+3);
        if((q = (quadnode *) hashtable_find(&quadtable, h)) != NULL) {
                                                        // check if pattern found in hash table is correct
            if( patt == q->nw &&  0ull == q->ne && 0ull == q->sw && 0ull ==q->se && !q->isnode) {
                q->activity++;q->lasttime=totsteps;
            }
            else {                                      // collision in hash table at leaf level
                npatt=patt*randomizer;
                h = patt_hash(npatt,npatt+1,npatt+2,npatt+3);
                if((q = (quadnode *) hashtable_find(&quadtable, h)) != NULL) {
                                                        // check if pattern found in hash table is correct
                    if( patt == q->nw &&  0ull == q->ne && 0ull == q->sw && 0ull ==q->se && !q->isnode) {
                        q->activity++;q->lasttime=totsteps;
                    }
                    else {                                      // collision in hash table at 8-leaf level
                        quadcollisions++;
                        fprintf(stderr,"at %d quadhash 2ndary 8-pattern collision %llx hash %llx collides with %llx %llx %llx %llx\n",
                                totsteps,patt,h,q->nw,q->ne,q->sw,q->se);
                    }
                }
                else {                                           // new node or pattern, save in hash table
                    q = hash_patt8_store(h,patt);
                }
            }
        }
        else {                                           // new node or pattern, save in hash table
            q = hash_patt8_store(h,patt);
        }
        nw = patt & 0xffff;ne = (patt>>16) & 0xffff;sw = (patt>>32) & 0xffff; se = (patt>>48) & 0xffff;
        if(nw) hash_patt4_find(nw);if(ne) hash_patt4_find(ne);   // find or store 8x8 64-bit subpatterns, updating activities and lasttime
        if(sw) hash_patt4_find(sw);if(se) hash_patt4_find(se);   // store if new, otherwise update
        return(q);
}
//.......................................................................................................................................................
extern inline uint64_t node_hash(const uint64_t a, const uint64_t b, const uint64_t c, const uint64_t d) {
// now not used : could use instead of patt_hash for hash_node_find
   uint64_t r = (65537*(d)+257*(c)+17*(b)+5*(a));
   r += (r >> 11);
   return r ;
}
//.......................................................................................................................................................
extern inline quadnode * hash_patt16_store(const uint64_t h, const uint64_t nw, const uint64_t ne, const uint64_t sw, const uint64_t se) {
    int nr1,nr1s;
    quadinit.hashkey = h;
    quadinit.isnode=0;
    quadinit.nw=nw;
    quadinit.ne=ne;
    quadinit.sw=sw;
    quadinit.se=se;
    quadinit.size= 16;
    quadinit.firsttime=totsteps;
    quadinit.lasttime=totsteps;
    POPCOUNT64C(nw,nr1);nr1s=nr1;
    POPCOUNT64C(ne,nr1);nr1s+=nr1;
    POPCOUNT64C(sw,nr1);nr1s+=nr1;
    POPCOUNT64C(se,nr1);nr1s+=nr1;
    quadinit.pop1s =nr1s;
    hashtable_insert(&quadtable, h,(quadnode *) &quadinit);
    return (quadnode *) hashtable_find(&quadtable, h);
}
//.......................................................................................................................................................
extern inline quadnode * hash_patt16_find(const uint64_t nw, const uint64_t ne, const uint64_t sw, const uint64_t se) {
        quadnode *q,*q8;
        uint64_t h,nnw,nne,nsw,nse;
        const uint64_t randomizer = 11400714819323198549ull;

        h = patt_hash(nw,ne,sw,se);
        if((q = (quadnode *) hashtable_find(&quadtable, h)) != NULL) {
                                                        // check if pattern found in hash table is correct
            if( nw == q->nw &&  ne == q->ne && sw == q->sw && se ==q->se && !q->isnode) {
                q->activity++;q->lasttime=totsteps;
            }
            else {                                      // collision in hash table at leaf level
                // quadcollisions++;
                // fprintf(stderr,"at %d quadhash pattern collision %llx %llx %llx %llx hash %llx collides %llx %llx %llx %llx\n",
                //    totsteps,nw,ne,sw,se,h,q->nw.patt,q->ne.patt,q->sw.patt,q->se.patt);
                nnw=nw*randomizer;
                nne=ne*randomizer;
                nsw=sw*randomizer;
                nse=se*randomizer;
                h = patt_hash(nnw,nne,nsw,nse);
                if((q = (quadnode *) hashtable_find(&quadtable, h)) != NULL) {
                                                        // check if pattern found in hash table is correct
                    if( nw == q->nw &&  ne == q->ne && sw == q->sw && se ==q->se && !q->isnode) {
                        q->activity++;q->lasttime=totsteps;
                    }
                    else {                                      // collision in hash table at leaf level
                        quadcollisions++;
                        fprintf(stderr,"at %d quadhash 2ndary pattern collision %llx %llx %llx %llx hash %llx collides %llx %llx %llx %llx\n",
                                totsteps,nw,ne,sw,se,h,q->nw,q->ne,q->sw,q->se);
                    }
                }
                else {                                           // new node or pattern, save in hash table
                    q = hash_patt16_store(h,nw,ne,sw,se);
                }
            }
        }
        else {                                           // new node or pattern, save in hash table
            q = hash_patt16_store(h,nw,ne,sw,se);
        }
        if(nw) q8 = hash_patt8_find(nw);if(ne) q8 = hash_patt8_find(ne);   // find or store 8x8 64-bit subpatterns, updating activities and lasttime
        if(sw) q8 = hash_patt8_find(sw);if(se) q8 = hash_patt8_find(se);   // store if new, otherwise update
        return(q);
}
//.......................................................................................................................................................
extern inline quadnode * hash_node_store(uint64_t h, uint64_t nw, uint64_t ne, uint64_t sw, uint64_t se) {
    quadnode *q;
    quadinit.hashkey = h;
    quadinit.isnode=1;
    quadinit.nw=nw;quadinit.ne=ne;quadinit.sw=sw;quadinit.se=se;
    quadinit.firsttime=totsteps;quadinit.lasttime=totsteps;
    quadinit.pop1s=0;quadinit.size=0;                           // incremented below
    quadinit.activity=1;quadinit.topactivity=0;                 // incremented only in quadimage at top level
    
    if (nw) {
        if ((q = (quadnode *) hashtable_find(&quadtable, nw)) == NULL) fprintf(stderr,"Error nw node not found in hashtable.\n");
        else {quadinit.pop1s+=q->pop1s;quadinit.size=q->size<<1;}
    }
    if (ne) {
        if ((q = (quadnode *) hashtable_find(&quadtable, ne)) == NULL) fprintf(stderr,"Error ne node not found in hashtable.\n");
        else {quadinit.pop1s+=q->pop1s;quadinit.size=q->size<<1;}
    }
    if (sw) {
        if((q = (quadnode *) hashtable_find(&quadtable, sw)) == NULL) fprintf(stderr,"Error sw node not found in hashtable.\n");
        else {quadinit.pop1s+=q->pop1s;quadinit.size=q->size<<1;}
    }
    if (se) {
        if((q = (quadnode *) hashtable_find(&quadtable, se)) == NULL) fprintf(stderr,"Error se node not found in hashtable.\n");
        else {quadinit.pop1s+=q->pop1s;quadinit.size=q->size<<1;}
    }
    
    hashtable_insert(&quadtable, h,(quadnode *) &quadinit);
    return (quadnode *) hashtable_find(&quadtable, h);
}
//.......................................................................................................................................................
extern inline quadnode * hash_node_find(const uint64_t nw, const uint64_t ne, const uint64_t sw, const uint64_t se) {
                                                                // should only be called if arguments are hashtable keys, not bit patterns
        quadnode *q;
        uint64_t h,nnw,nne,nsw,nse;
    
        h = node_hash(nw,ne,sw,se);                             // alternatively use patt_hash : but this is optimized for 64 bit keys (unlikely to have 0 or low values)
        if((q = (quadnode *) hashtable_find(&quadtable, h)) != NULL) {
            if(nw == q->nw && ne == q->ne && sw == q->sw && se == q->se && q->isnode) { // node found in hash table
                q->activity++;q->lasttime=totsteps;
                // Still need to update activities of 4 subnodes if non-zero !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            }
            else {                                              // collision in hash table at node level
                nnw=nw * 11400714819323198549ull;
                nne= ne * 11400714819323198549ull;
                nsw=(uint64_t) sw * 11400714819323198549ull;
                nse=(uint64_t) se * 11400714819323198549ull;
                h = node_hash(nnw,nne,nsw,nse);                 // alternatively use patt_hash
                if((q = (quadnode *) hashtable_find(&quadtable, h)) != NULL) {
                                                                // check if pattern found in hash table is correct
                    if(nw == q->nw && ne == q->ne && sw == q->sw && se == q->se && q->isnode) { // node found in hash table
                        q->activity++;q->lasttime=totsteps;
                    }
                    else {                                      // collision in hash table at leaf level
                        quadcollisions++;
                        fprintf(stderr,"at %d quadhash node 2ndary collision %llx %llx %llx %llx hash %llx collides %llx %llx %llx %llx\n", totsteps,
                                 nw, ne, sw, se, h, q->nw, q->ne,  q->sw, q->se);  // simple recording for now, later do further chaining or whatever
                    }
                }
                else q=hash_node_store(h,nw,ne,sw,se);          // new node, save in hash table
            }
        }
        else q=hash_node_store(h,nw,ne,sw,se);                  // new node, save in hash table

        return(q);
}
//.......................................................................................................................................................
uint64_t quadimage(uint64_t gol[], short unsigned int *patt, int log2n) {           // routine to generate a quadtree for an entire binary image of long words
                                                                                    // makes use of global linear and quadratic size variable N2 for sizing internal arrays golp,golq
                                                                                    // assumes that n is power of 2
    unsigned int ij,ij1,n3;
    uint64_t golp[N2>>6];
    quadnode * golq[N2>>8];
    int n = 1 << log2n;
    extern inline short unsigned int pack16neighbors(uint64_t gol[],int log2n);
    extern inline void pack64neighbors(uint64_t gol[],uint64_t golp[],int log2n);
    
    if(n<16) {                                                                      // n < 16
        if (n==8) {                                                                 // n == 8
            pack64neighbors(gol,golp,log2n);
            golq[0]=hash_patt8_find(golp[0]);
            golq[0]->topactivity++;
            return(golq[0]->hashkey);
        }
        else {                                                                      // n == 4,2,1   use smallpatts array to store patterns (more efficient than continued quadtree)
            *patt=pack16neighbors(gol,log2n);
            hash_patt4_find(*patt);
            smallpatts[*patt].topactivity++;
            return 0ull;                                                            // in this case image key is returned in *patt rather than quad key
        }
    }
    else  {                                                                         // n >= 16
        pack64neighbors(gol,golp,log2n);                                            // 8x8 blocks of gol pixels packed into single 64bit words in golp
        n3=n>>3;                                                                    // n3=n/8 is number of such 8x8 blocks along each side of square : n3 is at least 2 here (n>=16)
        // for(ij=0;ij<n3*n3;ij++) { if (ij%8 == 0) fprintf(stderr,"\n step %d ij %d",totsteps,ij);fprintf(stderr," %llx ",golp[ij]);} fprintf(stderr,"\n");
        quadtable.expansion_frozen = 1;                                             // freeze quad hash table against expansion ( to ensure valid pointers during array ops)
        for (ij=ij1=0;ij<n3*n3;ij+=2,ij1++) {                                       //  hash all 16x16 patterns (2x2 of golp words) found as leaves of the quadtree
            golq[ij1]=hash_patt16_find(golp[ij],golp[ij+1],golp[(ij+n3)],golp[(ij+n3)+1]); // hash_patt16_find(nw,ne,sw,se) adds quad leaf entry if pattern not found
            ij+= ((ij+2)&(n3-1)) ? 0 : n3;                                          // skip odd rows since these are the northern parts of quads generated on even rows
        }

        for(n3 >>= 1; n3>1; n3 >>= 1) {                                             // proceed up the hierarchy amalgamating 2x2 blocks to next level until reach top
            for (ij=ij1=0;ij<n3*n3;ij+=2,ij1++) {                                   // hash_node_find(nw,ne,sw,se) adds quad node entry if node not found & updates activities
                golq[ij1]=hash_node_find(golq[ij]->hashkey,golq[ij+1]->hashkey,golq[ (ij+n3)]->hashkey,golq[(ij+n3)+1]->hashkey);
                ij+= ((ij+2)&(n3-1)) ? 0 : n3;                                      // skip odd rows since these are the southern parts of quads generated on even rows
            }
        }
        // if(golq[0]!=NULL) if(golq[0]->activity > 1) fprintf(stderr,"step %d image already found at t = %d activity %d\n",totsteps,golq[0]->firsttime,golq[0]->activity);
        golq[0]->topactivity++;
        quadtable.expansion_frozen = 0;                                                         // unfreeze quad hash table
        return(golq[0]->hashkey);
    }
}
//.......................................................................................................................................................
int labelimage(uint64_t hashkeypatt, short unsigned int labelimg[], short unsigned int label, int offset) { // rebuild image from quadimage at with label
    short unsigned int patt;
    int n;
    quadnode * q;
    // uint64_t nw,ne,sw,se;
    extern inline void unpack16neighbors(const short unsigned golpw, short unsigned int labelimg[],const unsigned int label,const int offset);
    extern inline void unpack64neighbors(const uint64_t golpw, short unsigned int labelimg[],const unsigned int label,const int offset);
    if(hashkeypatt<65536) {patt = hashkeypatt; unpack16neighbors(patt,labelimg,label,offset);}
    else if((q = (quadnode *) hashtable_find(&quadtable, hashkeypatt)) != NULL) {
        if(q->isnode) {                                     // hashed item is regular quadnode
            n = q->size >> 1;
            if (q->nw) labelimage(q->nw, labelimg, label, offset);
            if (q->ne) labelimage(q->ne, labelimg, label, offset-(offset&Nmask)+((offset+n)&Nmask));
            if (q->sw) labelimage(q->sw, labelimg, label, (offset+n*N)&N2mask);
            if (q->se) labelimage(q->se, labelimg, label, (offset-(offset&Nmask)+((offset+n)&Nmask)+n*N)&N2mask);
        }
        else {                                              // hashed item is pattern
            n=8;
            if (q->nw) {if (q->nw <65536) {patt = q->nw; unpack16neighbors(patt,labelimg,label,offset);}
                        else unpack64neighbors(q->nw,labelimg,label,offset);}
            if (q->ne) {if (q->ne< 65536) {patt = q->ne; unpack16neighbors(patt,labelimg,label,offset-(offset&Nmask)+((offset+n)&Nmask));}
                        else unpack64neighbors(q->ne,labelimg,label,offset-(offset&Nmask)+((offset+n)&Nmask));}
            if (q->sw) {if (q->sw< 65536) {patt = q->sw; unpack16neighbors(patt,labelimg,label,(offset+n*N)&N2mask);}
                        else unpack64neighbors(q->sw,labelimg,label,(offset+n*N)&N2mask);}
            if (q->se) {if (q->se< 65536) {patt = q->se; unpack16neighbors(patt,labelimg,label,(offset-(offset&Nmask)+((offset+n)&Nmask)+n*N)&N2mask);}
                        else unpack64neighbors(q->se,labelimg,label,(offset-(offset&Nmask)+((offset+n)&Nmask)+n*N)&N2mask);}
        }
    }
    return label;
}
//----------------------------------------------------------- pack012,3neighbors ------------------------------------------------------------------------
#define deltaxy(ij,x,y)  ((ij - (ij&Nmask) + (((ij+(x))&Nmask) + (y)*N)) & N2mask)
#define deltaxyn(ij,x,y,log2n)  ((ij - (ij&((1<<log2n)-1)) + (((ij+(x))&Nmask) + (y)*N)) & N2mask)
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
extern inline void pack49neighbors(uint64_t gol[],uint64_t golp[], int nbhood) {              // routine to pack all up to 3rd neighbours in single word
    unsigned int ij,k;
    int nbx[6] = {1,0,2,0,-4,-4};
    int nby[6] = {0,1,0,2,0,-4};
    /* 48 49 52 53 32 33 36 37                                                   // indexing order of bits in packed 64-bit word for neighborhoods around 0 bit
       50 51 54 55 34 35 38 39                                                   // unused bits are masked out after six step recursive packing
       56 57 60 61 40 41 44 45
       58 59 62 63 42 43 46 47
       16 17 20 21  0  1  4  5
       18 19 22 23  2  3  6  7
       24 25 28 29  8  9 12 13
       26 27 30 31 10 11 14 15 */

    for (ij=0;ij<N2;ij++) golp[ij] = gol[ij];                                     // copy 1 bit gol to golp
    for(k=0;k<6;k++)                                                              // hierarchical bit copy and swap
        for (ij=0;ij<N2;ij++)
             golp[ij] |= golp[deltaxy(ij,nbx[k],nby[k])]<<(1<<k);                 // 8x8 packed arrays
    if(nbhood == 7)                                                               // masks out 15 values in top row and left column to give 7x7 neighbourhoods
        for (ij=0;ij<N2;ij++) golp[ij] = golp[ij]&0xfac8ffccfafaffffull;          // mask removes bit numbers 16,18,24,26,32,33,36,37,48,49,50,52,53,56,58
    else if (nbhood == 5)                                                         // masks in 25 values to give 5x5 neighbourhoods
        for (ij=0;ij<N2;ij++) golp[ij] = golp[ij]&0xf0005f0030f0135full;          // mask in bit nrs 0,1,2,3,4,6,8,9,12,20,21,22,23,28,29,40,41,42,43,44,46,60,61,62,63
    else                                                                          // masks in 9 values to give 3x3 neighbourhoods
        for (ij=0;ij<N2;ij++) golp[ij] = golp[ij]&0x80000c0000a0000full;          // mask in bit nrs 0,1,2,3,21,23,42,43,63

}
//.......................................................................................................................................................
extern inline short unsigned int pack16neighbors(uint64_t wgol[], int log2n) {    // routine to pack up to 4x4 subarray of binary square array wgol (nxn) into single word
    uint64_t golp;                                                      // assuming side length of square is n (power of 2)
    if (log2n==0) return((short unsigned int) wgol[0]);
    else if (log2n==1) return((short unsigned int) (wgol[0]+(wgol[1]<<1)+(wgol[2]<<2)+(wgol[3]<<3)));
    else if (log2n==2) {
        golp=0ull;
        golp |=  (wgol[0]+(wgol[1]<<1)+(wgol[4]<<2)+(wgol[5]<<3));
        golp |= ((wgol[2]+(wgol[3]<<1)+(wgol[6]<<2)+(wgol[7]<<3))<<4);
        golp |= ((wgol[8]+(wgol[9]<<1)+(wgol[12]<<2)+(wgol[13]<<3))<<8);
        golp |= ((wgol[10]+(wgol[11]<<1)+(wgol[14]<<2)+(wgol[15]<<3))<<12);
        return((short unsigned int) golp);
    }
    else {fprintf(stderr,"pack16neighbours called with not permitted value of log2n %d\n",log2n);return(0);}
}
//.......................................................................................................................................................
extern inline void unpack16neighbors(const short unsigned golpw, short unsigned int labelimg[],const unsigned int label,const int offset){
    int k,ij;
    if (golpw < 2) labelimg[0] = golpw ? label : 0;
    else if (golpw < 16) {
        for(k=0;k<4;k++) {
            ij = deltaxy(offset,k&0x1,k>>1);
            labelimg[ij] = (golpw>>k)&0x1 ? label : 0;
        }
    }
    else {
        for(k=0;k<16;k++) {
            int k1 = k&0x1; int k2 = (k>>1)&0x1;
            ij = deltaxy(offset,k1+(((k>>2)&0x1)<<1),k2+((k>>3)<<1));
            labelimg[ij] = (golpw>>k)&0x1 ? label : 0;
        }
    }
}
//.......................................................................................................................................................
extern inline int log2size(const short unsigned int golpw) {
    if (golpw < 2) return 0;
    else if (golpw < 16) return 1;
    else return 2;
}
//.......................................................................................................................................................
extern inline void pack64neighbors(uint64_t gol[],uint64_t golp[],int log2n) {    // routine to pack 8x8 subarrays of full binary array gol into single words
    int n = 1 << log2n;
    int ij,ij1,k;                                                                 // assuming golp length >= (n*n)>>2^6
    int n2 = n*n;

    if (log2n<3) {
        fprintf(stderr,"Error trying to pack too small an array into 64 bit words: need >= 8x8 have %d x %d\n",n,n);
        golp[0]=0;
        return;
    }
    for (ij1=0;ij1<(n2>>6);ij1++) golp[ij1]=0ull;
    for(k=0;k<64;k++) {
        int kx = ((k&0xf)&3)+(((k>>4)&1)<<2); int ky = ((k&0xf)>>2)+(((k>>4)>>1)<<2);
        for (ij1=0,ij=0;ij<n2;ij+=8) {
            // golp[ij1++] |= gol[(ij &(n-1))+(k&0x7)+n*((ij>>log2n)+(k>>3))]<<k;// blocked as 8*8 not compatible with smallpatts
            golp[ij1++] |= gol[(ij &(n-1))+ kx +n*((ij>>log2n)+ky)]<<k;          // blocked as 4*4*4
            ij+= ((ij+8)&(n-1)) ? 0 : (8-1)*n;
        }
    }
}
//.......................................................................................................................................................
extern inline void unpack64neighbors(const uint64_t golpw, short unsigned int labelimg[],const unsigned int label,const int offset){
    int k,ij;                                                                   // only unpacks one word
    for(k=0;k<64;k++) {                                                         // bits blocked as 4*4*4 so that 1st 16 bits are nw 4*4 block
        int k1 = k&0xf; int k2 = k>>4;
        ij = deltaxy(offset,(k1&3)+((k2&1)<<2),(k1>>2)+((k2>>1)<<2));
        labelimg[ij] = (golpw>>k)&0x1 ? label : 0;
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
    
    d=0;
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
    pack49neighbors(newgol,working,it_nbhood);                      // 3x3,5x5 or 7x7 packed newgol values in working
    if(offdx==0 && offdy==0 && offdt==0) {
        pack49neighbors(gol,golmix,it_nbhood);                      // 3x3,5x5 or 7x7 packed gol values in golmix
        compare_all_neighbors(golmix,working);                      // compare all 8 directions N E S W NE SE SW NW;
    }                                                               // output=golmix will contain packed numbers of 7x7 differences for all 8 directions
    else {
        if (offdt<=-maxPlane) offdt=-maxPlane;
        if(offdt>0) offdt = 0;
        pack49neighbors(planesg[(newPlane-offdt)%maxPlane],golmix,it_nbhood);
        compare_neighbors(golmix,working,offdx,offdy);              // compare with a single direction (north) for gliders
    }
}
//----------------------------------------------------------- connected component labelling -------------------------------------------------------------
extern inline short unsigned int lab_union(equivrec eqv[], short unsigned int i, short unsigned int j) {
// Combine two trees containing node i and j. Union by rank with halving - see https://en.wikipedia.org/wiki/Disjoint-set_data_structure
// Return the root of union tree in equivalence relation's disjoint union-set forest
    short unsigned int root = i;
    while (eqv[root].pt<root) {
        eqv[root].pt = eqv[eqv[root].pt].pt;        // halving algorithm to compress path on the fly
        root = eqv[root].pt;
    }
    if (i != j) {
        short unsigned int rooti = root;
        short unsigned int rootj;
        root = j;
        while (eqv[root].pt<root) {
            eqv[root].pt = eqv[eqv[root].pt].pt;   // halving algorithm to compress paths on the fly
            root = eqv[root].pt;
        }
        if (root == rooti) return root;
        rootj=root;
        // if (eqv[rooti].rank < eqv[rootj].rank) {                 // swap rooti, rootj so that rooti has larger rank : disrupts ordering!
        if (rooti > rootj) {                 // swap rooti, rootj so that rooti is lower of two
            rootj=rooti;
            rooti=root;
        }
        eqv[rootj].pt = rooti;                                 // merge rootj into rooti
        if (eqv[rooti].rank == eqv[rootj].rank) {
            eqv[rooti].rank++;
        }
    }
    return root;
}
//.......................................................................................................................................................
extern inline short unsigned int label_cell_Wu(int ij,short unsigned int *nlabel) {  // first version, slightly less efficient?, not used
    short unsigned int clabel,clabel1,labelij;
    clabel=label[deltaxy(ij,0,-1)];                  // deltaxy takes periodic BCs into account (b in Suzuki Decision Tree, see Wu et.al.)
    if(clabel) labelij=eqv[clabel].pt;         // b
    else {                                           // N (=b) unlabelled
        clabel=label[deltaxy(ij,1,-1)];              // (c in Suzuki DT)
        if(clabel) {
            clabel1=label[deltaxy(ij,-1,-1)];        // NW (a in Suzuki DT)
            if(clabel1) {
                labelij=lab_union(eqv,clabel,clabel1); // c,a
            }
            else {
                clabel1=label[deltaxy(ij,-1,0)];     // W (d in Suzuki DT)
                if(clabel1) {
                    labelij=lab_union(eqv,clabel,clabel1); // c,d
                }
                else {
                    labelij=eqv[clabel].pt;    // c
                }
            }
        }
        else {
            clabel=label[deltaxy(ij,-1,-1)];         // NW (a in Suzuki DT)
            if(clabel) {
                labelij=eqv[clabel].pt;        // a
            }
            else {
                clabel=label[deltaxy(ij,-1,0)];      // W (d in Suzuki DT)
                if(clabel) {
                    labelij=eqv[clabel].pt;    // d
                }
                else {
                    clabel=++*nlabel;                // new label if a,b,c,d all without labels
                    eqv[clabel].pt=clabel;
                    labelij=clabel;
                    eqv[clabel].rank=0;
                }
            }
        }
    }
    return labelij;
}
//.......................................................................................................................................................
extern inline short unsigned int label_cell(int ij,short unsigned int *nlabel) {
    short unsigned int clabel,clabel1,labelij,labelijold;
    labelijold=label[ij];
    clabel=label[deltaxy(ij,0,-1)];                  // deltaxy takes periodic BCs into account (b in Suzuki Decision Tree, 2010)
    if(clabel) labelij=eqv[clabel].pt;         // copy b
    else {                                           // N (=b) unlabelled
        clabel=label[deltaxy(ij,-1,0)];              // W (d in Suzuki DT)
        if(clabel) {
            clabel1=label[deltaxy(ij,1,-1)];         // NE (c in Suzuki DT)
            if(clabel1) {
                labelij=lab_union(eqv,clabel1,clabel); // resolve c,d
            }
            else {
                labelij=eqv[clabel].pt;        // copy d
            }
        }
        else {
            clabel=label[deltaxy(ij,1,-1)];              // NE (c in Suzuki DT)
            if(clabel) {
                clabel1=label[deltaxy(ij,-1,-1)];        // NW (a in Suzuki DT)
                if(clabel1) {
                    labelij=lab_union(eqv,clabel,clabel1); // resolve c,a
                }
                else {
                    labelij=eqv[clabel].pt;        // copy c
                }
            }
            else {
                clabel1=label[deltaxy(ij,-1,-1)];        // NW (a in Suzuki DT)
                if(clabel1) {
                    labelij=eqv[clabel1].pt;       // copy a
                }
                else {
                    clabel=++*nlabel;                    // new label if a,b,c,d all without labels
                    eqv[clabel].pt=clabel;
                    labelij=clabel;
                    eqv[clabel].rank=0;
                }
            }
        }
    }
    if (labelijold) labelij= lab_union(eqv,labelij,labelijold); // resolve labelij,labelijold needed for periodic BCs wrap around
    return labelij;
}
//.......................................................................................................................................................
void checklabels(equivrec eqv[],unsigned short int *nlabel) {
    short unsigned int xlabel = 0;
    short unsigned int i;
    for (i = 1; i < *nlabel + 1; i++) {
        if (eqv[i].pt > i) {
            fprintf(stderr,"Error in label equivalencies at t=%d for i=%d, parent %d\n",totsteps,i,eqv[i].pt);
            xlabel++;
        }
    }
    if (xlabel) fprintf(stderr,"Error for %d cases\n",xlabel);
}
//.......................................................................................................................................................
void flattenlabels(equivrec eqv[],unsigned short int *nlabel) {
// Flatten the Union-Find tree and relabel the components.
    short unsigned int xlabel = 1;
    short unsigned int i;
    for (i = 1; i <= *nlabel; i++) {
        if (eqv[i].pt < i) {
            eqv[i].pt = eqv[eqv[i].pt].pt;
        }
        else {
            eqv[i].pt = xlabel++;
        }
    }
    *nlabel = xlabel-1;
}
//.......................................................................................................................................................
short unsigned int label_components(uint64_t gol[]) {
// fast component labelling, updating global equivalence table equiv and placing labels in global array label, return number of labels
    int i,ij,sum,sumoverlap,maxoverlap,overallmaxoverlap;   //   abc   scan mask to calculate label at position e
    short unsigned int nlabel = 0;                          //   de
    short unsigned int oldnlabel;                           //   with lapmod need also: ret
    short unsigned int lab,conn,connprev,connf,maxlabel;
    float rsumoverlap;
    int nunique,nzconnect,nremconnect;
    int dx[9]={0,-1,0,1,1,1,0,-1,-1};
    int dy[9]={0,-1,-1,-1,0,1,1,1,0};
    static int connectout = 0;                              // whether to print lists of connected component mappings t-1 t
    void testflow(void);
    extern int maxmatch(int m, short unsigned int kk[], unsigned int ii[], short unsigned int xlap[], short unsigned int ylap[], short unsigned int dist[]);
    static int first = 1;
    if(!first) {
        oldnlabel = oldncomponents = ncomponents;
        for(ij=0;ij<N2;ij++) oldlabel[ij]=label[ij];
        for(i=0;i<NLM;i++) oldrelabel[i]=0;
        for(i=0;i<ncomponents;i++) oldrelabel[i]=relabel[i];
    }
    else {
        for(ij=0;ij<N2;ij++) oldlabel[ij]=0;
        oldnlabel = oldncomponents = 0;
        for(i=0;i<NLM;i++) oldrelabel[i]=0;
        first = 0;
    }
    for(ij=0;ij<N2;ij++) label[ij]=0;
    for(i=0;i<(NLM);i++) eqv[i].pt=0;
    for(i=0;i<(NLM);i++) eqv[i].rank=0;
    for(i=0;i<(NLM);i++) eqv[i].size=0;
    for(i=0;i<(NLM);i++) xlap[i]=ylap[i]=0;

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
    // checklabels(eqv,&nlabel);                            // check labels point to lower values in equivalence table
    flattenlabels(eqv,&nlabel);                             // single pass flatten of equivalence table

    for(ij=0;ij<N2;ij++) {                                  // do second pass of main array assigning final root labels and calculating component sizes in pixels
        if (gol[ij]) {
            label[ij]=eqv[label[ij]].pt;
            eqv[label[ij]].size++;
       }
    }

    for(i=0;i<NLM;i++) connlists[i]=0;                      // intialize all connection lists to zero
    for(i=0;i<NLM;i++) connlen[i]=0;
    for(i=0;i<NLM;i++) connlistsf[i]=0;
    for(i=0;i<NLM;i++) connlenf[i]=0;
    for(ij=0;ij<N2;ij++) {                                  // these are the open memory reserve of connection elements to be linked, better to use memset perhaps
        connections[ij].next=connections[ij].nextf=0;
        connections[ij].oldlab=connections[ij].newlab=0;
        connections[ij].overlap=0;
    }
    connused=0;

    for(ij=0;ij<N2;ij++) {                                  // build up backwards connections for each label[ij] via nbs in t-1 frame, insert in increasing label nr
        if(label[ij]) {
            for(int k=0;k<9;k++) {
                if((lab=oldlabel[deltaxy(ij,dx[k],dy[k])])) {
                    conn=connlists[label[ij]];
                    connprev=0;
                    while(conn && connections[conn].oldlab<lab) {
                        connprev = conn;
                        conn=connections[conn].next;
                    }
                    if (conn) {                             // connections[conn].oldlab>=lab : if oldlab>lab insert & increment overlap, else just increment overlap
                        if(connections[conn].oldlab>lab) {  // insert new node
                            if(connused>=N2) fprintf(stderr,"Error, out of connection memory\n");
                            connections[connused].oldlab=lab;
                            connections[connused].newlab=label[ij];
                            connections[connused].next=conn;
                            connections[connused].overlap++;
                            if(connprev) connections[connprev].next = connused++;
                            else connlists[label[ij]]=connused++;
                        }
                        else connections[conn].overlap++;
                    }
                    else {                                  // insert as last node in connection list added
                        if(connused>=N2) fprintf(stderr,"Error, out of connection memory\n");
                        connections[connused].oldlab=lab;
                        connections[connused].newlab=label[ij];
                        connections[connused].overlap++;
                        // connections[connused].next=0;
                        if(connprev) connections[connprev].next = connused++;
                        else connlists[label[ij]]=connused++;
                    }
                }
            }
        }
    }

    for(i=1;i<=nlabel;i++) {                                 // count number of connections to old components for each component
        sum = 0;
        conn=connlists[i];
        while(conn) {
            sum++;
            conn=connections[conn].next;
        }
        connlen[i]=sum;
    }

    for(i=1;i<=nlabel;i++) {                                 // count number of connections from each old component to new components
        conn=connlists[i];
        while(conn) {
            connlenf[connections[conn].oldlab]++;
            conn=connections[conn].next;
        }
    }
    connpref[0]=0;
    for(i=1;i<=nlabel;i++) {                                 // weave forward connections from each old component to new components & record preference
        conn=connlists[i];
        maxlabel=0; maxoverlap = 0;
        while(conn) {
            lab=connections[conn].oldlab;
            connf = connlistsf[lab];
            connprev=0;
            while(connf && connections[connf].newlab<i) {   // change to < for compatibility with lapmod sparse array order
                connprev = connf;
                connf = connections[connf].nextf;
            }
            if (connf) {                                    // connections[connf].newlab<=i : do nothing if newlab==i, otherwise insert
                if(connections[connf].newlab>i) {           // insert node conn
                    connections[conn].nextf=connf;
                    if(connprev) connections[connprev].nextf = conn;
                    else connlistsf[lab]=conn;
                }
            }
            else {                                          // insert as last node in connection list added
                if(connprev) connections[connprev].nextf = conn;
                else connlistsf[lab]=conn;
            }
            if(connections[conn].overlap>maxoverlap) {
                maxoverlap=connections[conn].overlap;
                maxlabel = connections[conn].oldlab;
            }
            conn=connections[conn].next;
        }
        connpref[i]=maxlabel;
    }

    overallmaxoverlap = 0;
    connpreff[0]=0;
    for(i=1;i<=oldnlabel;i++) {                             // calclulate weighted overlaps and max overlap connection and overall max overlap
        sumoverlap = 0;
        connf=connlistsf[i];
        maxlabel=0; maxoverlap = 0;
        while(connf) {
            sumoverlap+=connections[connf].overlap;
            if (connections[connf].overlap>maxoverlap) {
                maxoverlap=connections[connf].overlap;
                maxlabel = connections[connf].newlab;
            }
            connf=connections[connf].nextf;
        }
        if(maxoverlap>overallmaxoverlap) overallmaxoverlap = maxoverlap;
        connpreff[i]=maxlabel;
        
        rsumoverlap = 1.0/(float) sumoverlap;
        connf=connlistsf[i];                                // save weighted forward overlaps : .woverlap
        while(connf) {
            connections[connf].woverlap=rsumoverlap*connections[connf].overlap;
            connf=connections[connf].nextf;
        }
    }
    
    for(i=1;i<=nlabel;i++) {                                 // construct alternative overlap aoverlap from woverlap
        rsumoverlap = 0.;                                    // reuse rsumoverlap as floating point sum before doing reciprocal
        conn=connlists[i];
        while(conn) {
            rsumoverlap+=connections[conn].woverlap;
            conn=connections[conn].next;
        }
        rsumoverlap = 1.0/rsumoverlap;
        conn=connlists[i];                                  // save alternative weighted backward overlaps : .aoverlap
        while(conn) {
            connections[conn].aoverlap=rsumoverlap*connections[conn].woverlap;
            conn=connections[conn].next;
        }
    }

    for(i=1;i<=nlabel;i++) {                               // prune other connections for mutually preferred matchings
        if(connpreff[connpref[i]]==i) {
            conn=connlists[i];
            while(conn) {
                connlenf[connections[conn].oldlab]++;
                conn=connections[conn].next;
            }
        }
    }

    for(i=0;i<N2;i++) {                                    // reset arrays used in LAP
        cclap[i]=0;
        kklap[i]=0;
    }
    for(i=0;i<NLM;i++) iilap[i]=0;
    nlap = 0;                                              // initialize counters
    nclap = 0;
    iilap[0] = 0;
    for(i=1;i<=oldnlabel;i++) {                            // setup cost matrix as sparse array for lapmod
        if(connpref[connpreff[i]]==i) {                    // use only this connection, pruning all others if mutually preferred
            kklap[nclap]=connpreff[i];                     // for maxmatch, labels for sparse cost matrix column index kk are in range from 1 to nlabel
            cclap[nclap]=1;                                // minimal but non zero cost
            nclap++;
        }
        else {
            connf=connlistsf[i];
            while(connf) {
                lab = connections[connf].newlab;
                if(connpreff[connpref[lab]]!=lab) {           // only connect to labels that are not preassigned by mutuality to another label
                    kklap[nclap]=connections[connf].newlab;   // for maxmatch labels for sparse cost matrix column index kk are newlab : and run from 1 to nlabel
                    cclap[nclap]=connections[connf].overlap ?  (1+overallmaxoverlap-connections[connf].overlap) :  100 * overallmaxoverlap;   // use if minimizing cost
                    // if (cclap[nclap]<0) fprintf(stderr,"error in cost matrix from genelife, negative value at %d\n",nclap);
                    nclap++;
                }
                connf=connections[connf].nextf;
            }
        }
        nlap++;                                           // end of row, possibility of zero entries in row
        iilap[nlap]=nclap;
    }

    nmatched=maxmatch(nlap,kklap,iilap,xlap,ylap,dist);
                                                          // optionally, print connections, preferred, matched, and list of possible
    if (connectout) {
        fprintf(stderr,"BACKWARD\n");
        for(i=1;i<=nlabel;i++) {                          // print backward connections
            fprintf(stderr,"step %5d: %3d bwd conn's for newlabel %4d prefers %4d assigned y%4d:",totsteps,connlen[i],i,connpref[i],ylap[i]);
            conn = connlists[i];
            while(conn) {
                fprintf(stderr," %4d(%3d)",connections[conn].oldlab,connections[conn].overlap);
                conn=connections[conn].next;
            }
            fprintf(stderr,"\n");
        }
        fprintf(stderr,"FORWARD\n");
        for(i=1;i<=oldnlabel;i++) {                       // print forward connections
            fprintf(stderr,"step %5d: %3d fwd conn's for oldlabel %4d prefers %4d assigned x%4d:",totsteps,connlenf[i],i,connpreff[i],xlap[i]);
            conn = connlistsf[i];
            while(conn) {
                fprintf(stderr," %4d(%3d)",connections[conn].newlab,connections[conn].overlap);
                conn=connections[conn].nextf;
            }
            fprintf(stderr,"\n");
        }
    }
    
    for(i=0;i<NLM;i++) relabel[i]=0;                    // now using matching to implement relabelling of new components
    if(totsteps==0) {
        fprintf(stderr,"totsteps 0 nlabel %d\n",nlabel);
        for(i=1;i<=nlabel;i++) relabel[i]=i;
    }
    else {
        for(i=1;i<=NLM;i++) queue_array[i] = 0;
        for(i=1;i<=nlabel;i++) if(ylap[i]) {
            relabel[i]=oldrelabel[ylap[i]];             // keep old label for these matched components
            queue_array[relabel[i]]=1;                  // mark this label as taken
        }
        for(ij=i=1;i<=nlabel;i++)  {
            if(!ylap[i]) {                              // if unmatched component
                while(queue_array[ij]) ij++;            // find next free label ij with relabel[ij]==0, i.e. not yet assigned
                relabel[i]=ij;
                ij++;
            }
        }
    }
    //  for(i=1;i<=nlabel;i++) fprintf(stderr,"relabel[%4d]=%4d ylap[%4d]=%4d\n",i,relabel[i],i,ylap[i]);


    for(nunique=0,i=1;i<=oldnlabel;i++) if (connpref[connpreff[i]] == i) nunique++;
    for(nzconnect=0,i=1;i<=oldnlabel;i++) if (!connlistsf[i]) nzconnect++;
    for(nremconnect=0,i=1;i<=oldnlabel;i++) if (iilap[i]==iilap[i-1]) nremconnect++;       // this includes those components with all connections removed by nunique pairs
    if(!(totsteps % 10) && colorupdate1 && connectedprints) {
        fprintf(stderr,"connected cpts:  %d(%d) matched(unique) & %d(%d) with no-residual(no) connections i.e. %d out of %d(%d) old(new) components\n",
                        nmatched,nunique,nremconnect,nzconnect,nmatched+nremconnect,oldnlabel,nlabel);
    }
    
    /* do this for lapmod not maxmatch : cclap needs to be of type cost_t not short unsigned int
    for(i=nlap;i<=nlabel;i++) {                         // if nlabel > oldnlabel, then fill out cost matrix further with dummy nodes
        cclap[nclap]=(cost_t) 100 * overallmaxoverlap;  // alternatively use LARGE (results in error);
        kklap[nclap++]=i;
        nlap++;
        iilap[nlap]=nclap;
    } */
    // testflow();  // maximum flow test case
    /*
    fprintf(stderr,"Linear Assignment Problem cc,kk output\n");
    for(i=0,ij=0;ij<nclap;ij++) {
        while(iilap[i]<=ij) {
            i++;
            fprintf(stderr,"\n");
        }
        fprintf(stderr,"ij %d i %d j=kk[ij] %d cc[ij] %f ii[i] %d connlenf[i] %d\n",ij,i,kklap[ij],cclap[ij],iilap[i],connlenf[i]);
    }
    */


//    fprintf(stderr,"step %d lapmod called with nlap %d (out of %d) and nclap %d and ncomponents %d and maxvalue %d\n",totsteps,nlap,oldnlabel,nclap,nlabel,overallmaxoverlap);
//    ret=lapmod_internal(nlap,cclap,iilap,kklap,xlap,ylap,FP_DYNAMIC);
//    fprintf(stderr,"step %d lapmod returns %d with first three entries %d %d %d\n",totsteps,ret,xlap[0],xlap[1],xlap[2]);

    return nlabel;
}
//.......................................................................................................................................................
short unsigned int extract_components(uint64_t gol[]) {
    int i,j,ij,ij1,log2n;
    short unsigned int k,nside,nlabel,patt;
    short unsigned int histside[log2N+1];
    // uint64_t fullpatt;   // debugging
    int wpixels;
    uint64_t hashkey,rand;
    
    for(i=1;i<=ncomponents;i++) {
        oldcomplist[i]=complist[i];                                                                     // copy whole component struct to previous time step list
    }
    oldncomponents=ncomponents;
    ncomponents = label_components(gol);                                                                // label connected components at this time step
    nlabel = ncomponents;

    for(i=1;i<=nlabel;i++) {                                                                            // initialize component structures for horizontal scan
        complist[i].lastrc=0;
        complist[i].label=0;
        complist[i].pixels=0;
    }
    for (i=0;i<N;i++) {                                                                                 // find lateral limits of each component in horizontal scan
      for (j=0; j<N; j++) {
        ij = j*N+i;
        if(label[ij]>nlabel) fprintf(stderr,"in extract_components step %d label %d out of bounds\n",totsteps,label[ij]);
        if (label[ij]) {                                                                                // if site labelled
            if (!complist[label[ij]].label) {                                                           // if label encountered for first time
                complist[label[ij]].label=label[ij];
                complist[label[ij]].E=complist[label[ij]].W=i;                                          // set both horizontal bounds to first encountered x for label
                complist[label[ij]].lastrc=i;                                                           // row contains this label
            }
            else if (complist[label[ij]].lastrc != i) {                                                 // label reencountered for first time in row
                if (((complist[label[ij]].lastrc+1)&(N-1)) == i) {                                      // continuation of component from previous row
                    if(complist[label[ij]].E>=complist[label[ij]].W) complist[label[ij]].E=i;
                }
                else if (((complist[label[ij]].lastrc+1)&(N-1)) <  i) {                                 // component resumes after row gap
                    // fprintf(stderr,"HORIZ TRACK step %d setting W %d after gap at i %d j %d label %d lastrc %d\n",totsteps,i,i,j,label[ij],complist[label[ij]].lastrc);
                    complist[label[ij]].W=i;
                }
                complist[label[ij]].lastrc=i;                                                           // row contains this label
            }
            complist[label[ij]].pixels++;
        }
      }
    }
    for(i=1;i<=nlabel;i++) {                                                                            // initialize component structures for vertical scan
        complist[i].lastrc=0;
        complist[i].label=0;
    }
    for (j=0;j<N;j++) {                                                                                 // find vertical limits of each component
      for (i=0; i<N; i++) {
        ij = j*N+i;
        if(label[ij]>nlabel) fprintf(stderr,"in extract_components step %d label %d out of bounds\n",totsteps,label[ij]);
        if (label[ij]) {                                                                                // if site labelled
            if (!complist[label[ij]].label) {                                                           // if label encountered for first time
                complist[label[ij]].label=label[ij];
                complist[label[ij]].N=complist[label[ij]].S=j;                                          // set both vertical bounds to first encountered x for label
                complist[label[ij]].lastrc=j;                                                           // col contains this label
            }
            else if (complist[label[ij]].lastrc != j) {                                                 // label reencountered for first time in col
                if (((complist[label[ij]].lastrc+1)&(N-1)) == j) {                                      // continuation of component from previous col
                    if(complist[label[ij]].S>=complist[label[ij]].N) complist[label[ij]].S=j;
                }
                else if (((complist[label[ij]].lastrc+1)&(N-1)) <  j) {                                 // component resumes after col gap
                    // fprintf(stderr,"VERT  TRACK step %d setting N %d after gap at i %d j %d label %d lastrc %d\n",totsteps,j,i,j,label[ij],complist[label[ij]].lastrc);
                    complist[label[ij]].N=j;
                }
                complist[label[ij]].lastrc=j;                                                           // col contains this label
            }
        }
      }
    }

    for(i=1;i<=nlabel;i++) {
        nside=((complist[i].E-complist[i].W)&(N-1));                                                    // modulo calculation required if component crosses periodic boundary
        k = ((complist[i].S-complist[i].N)&(N-1));
        nside = nside > k ? nside+1 : k+1;                                                              // side length of square one more than greater of vertical & horiz. difference
        for(log2n=0;(1<<log2n)<nside;log2n++);
        complist[i].log2n = log2n;                                                                      // log2n is smallest power of 2 for side of square enclosing component
        nside = 1<<log2n;                                                                               // convert nside to next power of 2 for quadtree analysis;
        wpixels = 0;
        for(ij=0;ij<nside*nside;ij++) {
            ij1 =   ((complist[i].W+(ij&(nside-1)))&(N-1)) +                                            // calculate coordinate ij1 in full array of ij index in the component square
                 N*((complist[i].N+(ij>>log2n))&(N-1));
            working[ij] = (label[ij1]==i) ? 0x1ull : 0ull;                                              // use working array to store binary component image
            if (working[ij]) wpixels++;
        }
        hashkey = quadimage(working,&patt,log2n);                                                       // quadtree hash code of component image: needs to work for all nside=2^n
        if(!hashkey) {                                                                                  // components either have (quad!=NULL and patt==0) or (quad==NULL and patt!=0)
            complist[i].patt=patt;
            complist[i].quad=0ull;
            /* DEBUGGING
            fullpatt = (uint64_t) patt;
            POPCOUNT64C(fullpatt, ij1);
            if (ij1 != complist[i].pixels) {
                fprintf(stderr,"step %d label %d ERROR in small pattern construction patt pixels : patt %x pixels %d fullpatt %llx n1s %d wpixels = %d nside %d\n",
                                                    totsteps, i, patt, complist[i].pixels, fullpatt, ij1, wpixels, nside);
                fprintf(stderr,"                 N S W E = (%4d %4d %4d %4d) lastrc %d\n",complist[i].N,complist[i].S,complist[i].W,complist[i].E,complist[i].lastrc);
            }
            */
        }
        else {
            complist[i].patt=0;
            complist[i].quad=hashkey;
        }
    }

    
    if(totsteps==0) {                                                                                         // initialize colors to labels for t=0
        for(i=1;i<=nlabel;i++) {
            complist[i].gcolor=(float)(i-1)/ (float) nlabel;
        }
    }
    else {                                                                                              // for t>0 mixcolors from overlapping connected components
        for(i=1;i<=nlabel;i++) {
            RAND128P(rand);
            complist[i].gcolor=mixcolor(i,rand);
        }
    }
    
    for(i=0;i<=log2N;i++) {
        histside[i]=0;
    }
    for(i=1;i<=nlabel;i++) {
        histside[complist[i].log2n]++;
        /* fprintf(stderr,"step %d comp %d pixels %d rect (%d,%d, %d,%d)\n",totsteps,i,complist[i].pixels,complist[i].W,complist[i].N,complist[i].E,complist[i].S);
           if(complist[i].log2n>=8) {
            fprintf(stderr,"step %d component %d with log2n %d N,S,W,E %d %d %d %d\n",
                    totsteps,i,complist[i].log2n,complist[i].N,complist[i].S,complist[i].W,complist[i].E);
            for(ij=0;ij<N2;ij++) if(label[ij]==i) label[ij]=0xffff;                                    // recolour component white for debugging
        } */
    }

    // fprintf(stderr,"step %d histogram of component log2n\n",totsteps);
    if(!(totsteps % 10) && colorupdate1 && connectedprints) {
        fprintf(stderr,"histogram log2n ");for(i=0;i<=log2N;i++) fprintf(stderr," %5d",i);fprintf(stderr,"\n");
        fprintf(stderr,"conn cmpt counts");for(i=0;i<=log2N;i++) fprintf(stderr," %5d",histside[i]);fprintf(stderr,"\n");
    }

    return nlabel;
}
//----------------------------------------------------------------- geography ---------------------------------------------------------------------------
void v_scroll(uint64_t newgol[],uint64_t newgolg[],uint64_t newgolb[],uint64_t newgolr[]) {
    int ij,scroll_needed;

    scroll_needed = 0;
    for (ij=N2-N;ij<N2;ij++) {                              // clear top row only if needed
        if(newgol[ij]) {
            scroll_needed = 1;
            break;
        }
    }
    last_scrolled = scroll_needed;                          // global variable needed for correct glider detection if scrolling

    for (ij=0;ij<N;ij++) {                                  // delete genes in bottom buffer row
        if(newgol[ij]) {
            newgol[ij]=0ull;
            if(diagnostics & diag_hash_genes)
                hashdeletegene(newgolg[ij],newgolb[ij],"error in v_scroll hashdeletegene call for step %d with gene %llx\n");
            newgolg[ij]=gene0;
            newgolb[ij]=0ull;
            newgolr[ij]=0ull;
        }
    }

    if(!scroll_needed) return;

    for (ij=0;ij<N2-N;ij++) {                               // scroll rows down 1 leaving top row intact
        newgol[ij]=newgol[ij+N];
        newgolg[ij]=newgolg[ij+N];
        newgolb[ij]=newgolb[ij+N];
        newgolr[ij]=newgolr[ij+N];
    }
    for (ij=0;ij<N;ij++) {                                  // delete all states and genes in new bottom buffer row
        if(newgol[ij]) {
            newgol[ij]=0ull;
            if(diagnostics & diag_hash_genes)
                hashdeletegene(newgolg[ij],newgolb[ij],"error in v_scroll hashdeletegene call for step %d with gene %llx\n");
            newgolg[ij]=gene0;
            newgolb[ij]=0ull;
            newgolr[ij]=0ull;
        }
    }
    for (ij=N2-N;ij<N2;ij++) {                              // clear top row
        if(newgol[ij]) {
            newgol[ij]=0ull;
            if(diagnostics & diag_hash_genes)
                hashdeletegene(newgolg[ij],newgolb[ij],"error in v_scroll hashdeletegene call for step %d with gene %llx\n");
            newgolg[ij]=gene0;
            newgolb[ij]=0ull;
            newgolr[ij]=0ull;
        }
    }
}
//.......................................................................................................................................................
void random_influx(uint64_t newgol[],uint64_t newgolg[],uint64_t newgolb[],uint64_t newgolr[]) {
    int Nf,i,j,ij,i0,j0,i1,j1,d;
    uint64_t randnr,mask,parentid;
    static unsigned int rmask = (1 << 15) - 1;
    unsigned int density;
    parentid = ((uint64_t) totsteps)<<32;
    
    if(rbackground) {                                           // homogeneous random background at rate rbackground/32768 per site per frame
        if(randominflux<3) {                                    // only create new genes as perturbation if randominflux<3
            if(randominflux < 2 || selection < 8) {             // random gene if randominflux<2 or selection < 8
                for(ij=0;ij<N2;ij++) {
                    if((rand()&rmask) < rbackground) {          // random event
                        if (!newgol[ij]) {                      // if not live, random genome
                            newgol[ij] = 1ull;
                            RAND128P(randnr);
                            newgolg[ij] = randnr;
                            if(diagnostics & diag_hash_genes) hashaddgene(ij,newgolg[ij],rootgene,newgolb+ij,(parentid+rootclone+ij),0x1ull);
                            newgolr[ij]=0ull;
                        }
                    }
                }
            }
            else if(randominflux==2) {                          // randominflux==2 and selection>=8
                for(ij=0;ij<N2;ij++) {
                    if((rand()&rmask) < rbackground) {
                        if (!newgol[ij]) {                      // if not live, fill with game of life genome
                            newgol[ij] = 1ull;
                            newgolg[ij] = genegol[selection-8]; // GoL gene for particular selection model coding
                            if(diagnostics & diag_hash_genes) hashaddgene(ij,newgolg[ij],rootgene,newgolb+ij,(parentid+rootclone+ij),0x1ull);
                            newgolr[ij]=0ull;
                        }
                    }
                }
            }
        }
        else if (randominflux==3) {                             // deletion perturbations only
            for(ij=0;ij<N2;ij++) {
                if((rand()&rmask) < rbackground) {              // random event
                    if (newgol[ij]) {                           // if live cell, delete gene
                        newgol[ij]=0ull;
                        if(diagnostics & diag_hash_genes) hashdeletegene(newgolg[ij],newgolb[ij],"error in randominflux=3 hashdeletegene call for step %d with gene %llx\n");
                        newgolg[ij]=gene0;
                        newgolb[ij]=0ull;
                        newgolr[ij]=0ull;
                    }
                }
            }
        }
        else if (randominflux==4) {                             // deletion perturbations s-dependent
            int s, se, k;
            int nb[8], ij, i, j, jp1, jm1, ip1, im1;
            uint64_t gols, nb1i;
            for(ij=0;ij<N2;ij++) {
                if((rand()&rmask) < rbackground) {              // random event
                    if (newgol[ij]) {                           // if live cell, delete gene
                        // compute s
                        i = ij & Nmask;  j = ij >> log2N;                                          // row & column
                        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                          // toroidal (j+1)*N and (j-1)*N
                        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                                // toroidal i+1, i-1
                        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;                  // new order of nbs
                        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
                        for (s=se=0,nb1i=0ull,k=0;k<8;k++) {                                       // packs non-zero nb indices in first up to 8*4 bits
                            gols=newgol[nb[k]];                                                    // whether neighbor is alive
                            s += gols;                                                             // s is number of live nbs
                            se += k&0x1&gols;                                                      // se is number of edge-centred live neighbours (odd k)
                            nb1i = (nb1i << (gols<<2)) + (gols*k);                                 // nb1i is packed list of live neighbour indices
                        }
                        // compute s done
                        if(s>1){
                            newgol[ij]=0ull;
                            if(diagnostics & diag_hash_genes) hashdeletegene(newgolg[ij],newgolb[ij],"error in randominflux=4 hashdeletegene call for step %d with gene %llx\n");
                            newgolg[ij]=gene0;
                            newgolb[ij]=0ull;
                            newgolr[ij]=0ull;
                        }
                    }
                }
            }
        }
        return;
    }

    if(randominflux>=2)
        if ((totsteps & 0xf) || randominflux == 3) return;      // only execute remaining code once every 16 time steps and if randominflux!=3

    density = initial1density;
    mask = ~0ull;
    Nf = initfield;
    if (Nf==0 || Nf>N) Nf=N;
    i0 = j0 = (N>>1)-(Nf>>1);

    for (i=0; i<Nf; i++) {
        for (j=0; j<Nf; j++) {
            ij=i0+i+N*(j0+j);
            if(randominflux==2) {                                 // border feathering as well as intermittent every 16 steps
                i1 = i<j ? i : j;                               // swap so that i1<=j1
                j1 = i<j ? j : i;
                d= j1< (Nf>>1) ? i1 : (i1 < Nf-j1 ? i1 : Nf-j1);// find Manhatten distance to border ij1
                density = (d <= 8 ? (initial1density >> (16-(d<<1))) : initial1density);
            }
            if(!newgol[ij]) {
                newgol[ij] = ((rand() & rmask) < density)?1ull:0ull;
 
                if (newgol[ij]) {                               // if live cell, fill with game of life genome or random genome
                    if (selection<8) {
                        RAND128P(randnr);
                        newgolg[ij] = gene0^(randnr&mask);
                    }
                    else {
                        newgolg[ij] = genegol[selection-8];
                    }
                    if(diagnostics & diag_hash_genes) hashaddgene(ij,newgolg[ij],rootgene,newgolb+ij,(parentid+rootclone+ij),0x1ull);
                }
            }
        }
    }
}
//---------------------------------------------------------------- save and retrieve gol... data --------------------------------------------------------
int savegols( int step, uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[],uint64_t golr[]) {
    FILE *fp;
    char fname[30];

    sprintf(fname,"gol_gsbr_data%d.ext",step);
    fp = fopen( fname , "wb" );
    
    fwrite(gol  , sizeof(uint64_t) , N2 , fp );
    fwrite(golg , sizeof(uint64_t) , N2 , fp );
    fwrite(golgstats, sizeof(uint64_t) , N2 , fp );
    fwrite(golb , sizeof(uint64_t) , N2 , fp );
    fwrite(golr , sizeof(uint64_t) , N2 , fp );

    fclose(fp);
    return 0;
}
//---------------------------------------------------------------- save data ----------------------------------------------------------------------------
int retrievegols( int step, uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[],uint64_t golr[]) {
    FILE *fp;
    char fname[30];

    sprintf(fname,"gol_gsbr_data%d.ext",step);
    fp = fopen( fname , "rb" );
    
    fread(gol ,  sizeof(uint64_t) , N2, fp );
    fread(golg ,  sizeof(uint64_t) , N2, fp );
    fread(golgstats, sizeof(uint64_t) , N2, fp );
    fread(golb , sizeof(uint64_t) , N2, fp );
    fread(golr , sizeof(uint64_t) , N2, fp );

    fclose(fp);
    return 0;
}
//---------------------------------------------------------------- update_23 ----------------------------------------------------------------------------
void update_23(uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[],uint64_t golr[],uint64_t newgol[], uint64_t newgolg[], uint64_t newgolgstats[], uint64_t newgolb[],uint64_t newgolr[]){
    // update for dissection of genetic rule variants within nearest neighbor sum s=2 or 3 only for survival and birth
    // update GoL for toroidal field which has side length which is a binary power of 2
    // encode without if structures for optimal vector treatment
    int s, s0, k, k1, kodd, nmut, crot;
    unsigned int kch,rulemodij,add2nd,mask1st;
    unsigned int select23live,pos_canon_neutral,survival,overwrite,enforcebirth,add2ndmask1st,nongolnottwice;
    int nb[8], nbc, nbch, ij, i, j, jp1, jm1, ip1, im1, ij1;
    uint64_t g, gs, nb1i, nb2i, randnr, r2;
    uint64_t nbmask, nbmaskr;
    uint64_t newgene, ancestor, livegenes[3], parentid;
    uint64_t s2or3, birth, statflag, nextgolstate;
    // short unsigned int patt; // needed if quadimage uncommented

    canonical = repscheme & R_2_canonical_nb;
    survival = survivalmask;
    overwrite = overwritemask;
    select23live = ((repscheme & R_0_2sel_3live)?1:0)+((repscheme & R_1_2sel_2live)?2:0);
    pos_canon_neutral = ((repscheme & R_2_canonical_nb)?1:0)+((repscheme & R_3_neutral_pos)?2:0);
    enforcebirth =((repscheme & R_4_enforce3birth)?1:0)+((repscheme & R_5_enforce2birth)?2:0);
    add2ndmask1st = ((repscheme & R_6_2ndnb_genes)?1:0)+((repscheme & R_7_1stnb_masks)?2:0);
    nongolnottwice = ((repscheme & R_8_nongolstat)?1:0)+((repscheme & R_9_nongolstatnbs)?2:0);
    add2nd = add2ndmask1st&0x1;
    parentdies = (repscheme & R_14_parentdies_23) ? 1 : 0;

    if(parentdies) for (ij=0; ij<N2; ij++) newgolgstats[ij] = 0ull;         // need to update statistics of neighbours with parenting information, so init required

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
          if (repscheme & R_19_quadrant_2nb1)  add2ndmask1st =     (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
          if (repscheme & R_20_quadrant_ngol)  nongolnottwice =    (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
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
          statflag |= F_livenbs & (nbmask<<16);                             // record live neighbour pattern
        } // end if s>1
        if (s2or3) {                                                        // if 2 or 3 neighbours alive
            if (repscheme & R_quadrant) {                                   // quarter the plane with 4 different parameter values
                if (repscheme & R_16_quadrant_sele)  select23live =      (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
                if (repscheme & R_17_quadrant_posn)  pos_canon_neutral = (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
                if (repscheme & R_18_quadrant_enfb)  enforcebirth =      (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
                if (repscheme & R_21_quadrant_surv)  survival =          (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
                if (repscheme & R_22_quadrant_over)  overwrite =         (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
                canonical = pos_canon_neutral&0x1;       // global value since needed in ...difft2-6 subroutines
            }
            birth = 0ull;
            newgene = 0ull;
            parentid = 0ull;
            kch = 0;

            if (s&0x1ull) {  // s==3                                        // allow birth (with possible overwrite)
              statflag |= F_3_live;                                         // record instance of 3 live nbs
              if ((0x1ull&overwrite)||!gol[ij] ) {                          // central site empty or overwrite mode
                birth = 1ull;                                               // birth flag
                for(k=0;k<s;k++) livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]]; // live neighbour genes
                kch=selectdifft3(nbmask, &crot, &kodd);
                parentid=golb[nb[kch&0x7]];                                 // default parent is for kch unless selected below
                if((livegenes[0]^livegenes[1])|(livegenes[0]^livegenes[2])) { // genes not all same, need ancestor calculation
                  nbch=nb[kch];
                  if (select23live&0x1) {                                   // execute selective replication of one of two otherwise unchosen live genes
                      nb2i = 0ull;
                      for(k1=k=0;k<3;k++) {                                 // choice of two other live genes for possible ancestor
                          nbc=(nb1i>>(k<<2))&0x7;
                          if(nb[nbc]!=nbch) nb2i = (nbc<<(k1++<<2))+nb2i;
                      }
                      if (add2nd) selectone_nbs(s,nb2i,nb,gol,golg,golb,&birth,&newgene,&parentid,&kch);  //2nd nb modulation
                      else selectone_of_2(s,nb2i,nb,golg,golb,&birth,&newgene,&parentid,&kch);
                      if (birth==0ull) {                                    // optional reset of ancestor & birth if no ancestors chosen in selectone
                        if((enforcebirth&0x1)||rulemodij)  {                // birth cannot fail or genes don't matter or no modification to gol rules
                            newgene = golg[nbch];
                            parentid=golb[nbch];
                            birth = 1ull;
                        }
                      }
                      else statflag |= F_2select;                           // ancestor has been chosen in selectone_of_2
                  }
                  else {
                      newgene = golg[nbch];
                      parentid = golb[nbch];
                  }
                } // end if not all live neighbors the same
                else {
                    statflag |= F_3g_same;
                    newgene = livegenes[0];                                 // genes all the same : copy first one, parentid default
                    if((~enforcebirth&0x1) && rulemodij) birth = 0ull;      // no birth for 3 identical genes if not enforcebirth3 and rulemod
                    kch=selectdifft3(nbmask, &crot, &kodd);                 // need a single deterministic ancestor for genealogy and golr
                }
              } // end central site empty or overwrite mode
            }  // end if s==3
            else {  // s==2                                                 // possible birth as exception to GoL rule
                statflag |= F_2_live;
                if (((select23live>>1)&0x1)&&(rulemodij||gol[ij])) {        // rule departure from GOL allowed or possible overwrite
                    if ((0x1ull&(overwrite>>1))||!gol[ij]) {                // either overwrite on for s==2 or central site is empty
                        if (add2nd) {
                            selectone_nbs(s,nb1i,nb,gol,golg,golb,&birth,&newgene,&parentid,&kch); //2nd nb modulation
                        }
                        else {
                            nbmask = (0x1ull<<(nb1i&0x7)) + (0x1ull<<((nb1i>>4)&0x7));
                            kch=selectdifft2(nbmask, &crot, &kodd);
                            if ((pos_canon_neutral>>1)&0x1) {               // enforce gene independent birth for s = 2 (corrected 31.1.2019, remove repscheme&)
                                newgene = golg[nb[kch]];
                                parentid= golb[nb[kch]];
                                birth = 1ull;
                            }
                            else {
                                selectone_of_2(s,nb1i,nb,golg,golb,&birth,&newgene,&parentid,&kch);
                            }
                            if(repscheme & R_10_13_2birth_k4) {             // birth only on the active subset of 4 canonical 2-live nb configs if any are active
                                if(~repscheme & (R_10_2birth_k0<<(kch-1))) birth = 0ull;   // cancel birth if chosen configuration not active
                            }
                        }

                        if(!birth && (enforcebirth&0x2)) {
                            nbmask = (0x1ull<<(nb1i&0x7)) + (0x1ull<<((nb1i>>4)&0x7));
                            kch=selectdifft2(nbmask, &crot, &kodd);
                            newgene = golg[nb[kch]];
                            parentid= golb[nb[kch]];
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
                    if(r2) {
                        statflag = statflag | F_mutation;
                        statflag = statflag | F_survmut;
                    }
                }
                if(diagnostics & diag_hash_genes) {
                    if(gol[ij]) {                                               // central old gene present: overwritten
                        hashdeletegene(golg[ij],golb[ij],"step %d hash delete error 1 in update, gene %llx not stored\n");
                    }
                    hashaddgene(ij,newgene,ancestor,newgolb+ij,parentid,statflag & F_mutation);
                }
                newgol[ij]  =  1ull;                                        // new game of life cell value: alive
                newgolg[ij] =  newgene;                                     // if birth then newgene
                newgolr[ij] = (golr[nb[kch]]<<4) | kch | 0x8;               // register ancestor offset index in record_of_dynamics_gene golr along with 0x8 for birth event
                statflag = statflag | F_birth;
                if (parentdies) {
                    ij1 = nb[kch];
                    newgolgstats[ij1] = newgolgstats[ij1] | F_parent;
                }
            } // end birth
            else {
                if ((survival&s&0x1ull)|((survival>>1)&(~s)&0x1ull)|((~rulemodij)&0x1ull)) {// (surv bit 0 and s==3) or (surv bit 1 and s==2) or not rulemod1ij
                // if ((survival&s&0x1ull)|((survival>>1)&(~s)&0x1ull)) {   // survival bit 0 and s==3, or (survival bit 1 and s==2)
                    newgol[ij]  = gol[ij];                                  // new game of life cell value same as old
                    newgolg[ij] = golg[ij];                                 // gene stays as before, live or not
                    newgolb[ij] = golb[ij];
                    newgolr[ij] = (golr[ij]<<4) | (s-1);                    // register s-1 in record_of_dynamics_gene golr along with 0 value of bit3 for survival event
                    if(gol[ij]) {
                        statflag |= F_survival;
                        if (golgstats[ij]&F_survmut) statflag |= F_survmut;
                    }
                }
                else {
                    if(diagnostics & diag_hash_genes) {
                        if(gol[ij]) {                                       // death : need to update hash table
                            hashdeletegene(golg[ij],golb[ij],"step %d hash delete error 2 in update, gene %llx not stored\n");
                        }
                    }
                    newgol[ij]  = 0ull;                                     // new game of life cell value dead
                    newgolg[ij] = 0ull;                                     // gene dies or stays dead
                    newgolb[ij] = 0ull;                                     // clone removed from site
                    newgolr[ij] = 0ull;                                     // displacement history erased on death
                    if(gol[ij]) statflag |= F_death;
                }
            } // end no birth
        }  // end if s2or3
        else {                                                              // else not birth or survival, 0 values for gol and gene
            if(diagnostics & diag_hash_genes) {
                if(gol[ij]) {                                               // death : need to update hash table
                    hashdeletegene(golg[ij],golb[ij],"step %d hash delete error 3 in update, gene %llx not stored\n");
                }
            }
            newgol[ij]  = 0ull;                                             // new game of life cell value
            newgolg[ij] = 0ull;                                             // gene dies
            newgolb[ij] = 0ull;                                             // clone removed from site
            newgolr[ij] = 0ull;
            if(gol[ij]) statflag |= F_death;
        }
        if(gol[ij]) statflag |= F_golstate;
        if(newgol[ij]^nextgolstate) statflag |= F_notgolrul;
        if(gol[ij]^newgol[ij]) {
            statflag |= F_golchange;
            if(statflag&F_notgolrul) statflag |= F_nongolchg;
        }
        else if (golgstats[ij]&F_nongolchg) statflag |= F_nongolchg;        // maintain non-GoL chg status until state changed by GoL rule
        if (parentdies) newgolgstats[ij] |= statflag;                       // newgolgstats may already contain updated parenthood info F_parent
        else newgolgstats[ij] = statflag;
    }  // end for ij

    if(parentdies) {
        for (ij=0; ij<N2; ij++) {
            statflag = newgolgstats[ij];
            if(gol[ij] && !(statflag&F_birth) && !(statflag&F_death) && (statflag&F_parent)) {
                newgol[ij]  = 0ull;                                         // new game of life cell value dead
                newgolg[ij] = 0ull;                                         // gene dies
                newgolb[ij] = 0ull;
                newgolr[ij] = 0ull;
                if(diagnostics & diag_hash_genes)
                    hashdeletegene(golg[ij],golb[ij],"step %d hash delete error 2 in update, gene %llx not stored\n");
                newgolgstats[ij] |= F_parentaldeath;
            }
        }
    }

    if(randominflux) random_influx(newgol,newgolg,newgolb,newgolr);
    if(vscrolling) v_scroll(newgol,newgolg,newgolb,newgolr);
    if ((colorfunction == 8) || (colorfunction2 == 8)) packandcompare(newgol,working,golmix);
    if(diagnostics & diag_component_labels) ncomponents=extract_components(newgol);

    for (ij=0; ij<N2; ij++) {       // complete missing hash table records of extinction and activities
        if(gol[ij]) hashgeneextinction(golg[ij],"hash extinction storage error %d in update at step %d, gene %llx not stored\n");
        if(newgol[ij]) hashgeneactivity(newgolg[ij],"hash activity storage error in update, gene %llx not stored\n");
        if(newgolb[ij]) hashcloneactivity(newgolb[ij],"hash activity storage error in update, clone %llx not stored\n");
    }
    // if(diagnostics & diag_hash_patterns) qimage = quadimage(newgol,&patt,log2N); // quadtree hash of entire image
}

//---------------------------------------------------------------- finish_update ------------------------------------------------------------------------
// two inline subroutines that are in common for all LUT based update rules
extern inline void finish_update_ij(int ij,int s,uint64_t golij,uint64_t gols,uint64_t nb1i,uint64_t nbmask,int nb[],uint64_t survive,uint64_t birth,
                                    uint64_t gol[],uint64_t golg[],uint64_t golgstats[],uint64_t golb[],uint64_t golr[],
                                    uint64_t newgol[],uint64_t newgolg[],uint64_t newgolgstats[],uint64_t newgolb[],uint64_t newgolr[]) {
        uint64_t randnr,r2,newgene,parentid,ancestor,statflag;
        int k, nmut, kodd, crot, nbest, nsame;
        unsigned int kch;
        statflag = F_livenbs & (nbmask<<16);                                    // record live neighbour pattern
        statflag |= F_s_live & (s&0x7);                                         // requires F_s_live to be in lowest 3 bits : val 8 mapped to 0
        if(birth) {                                                         // birth allowed by rules encoded in local genes (may still fail by selection)
            if (repscheme & R_7_random_resln) {
                RAND128P(randnr);                                           // inline exp so compiler recognizes auto-vec,
                kch = ((randnr>>32)&0xffff) % s;                            // choose random in this option only
                ancestor = newgene = golg[(nb[(nb1i>>(kch<<2))&0x7])];
                parentid = golb[(nb[(nb1i>>(kch<<2))&0x7])];
            }
            else {
                if ((ancselectmask>>(s-1)) &0x1) {                          // use genes to select ancestor
                    nbest=selectone_of_s(&kch,s,nb1i,nb,golg,golb,golr,&birth,&newgene,&parentid,&nbmask,ij);// selection scheme depends on repscheme parameter, selection depends on genes
                    ancestor = newgene;
                    // if(ij==IJDEBUG) fprintf(stderr,"DEBUG Difference 1 at totsteps=%d ij=%d s=%d nbest=%d birth=%llx newgene=%llx nb1i=%llx kch=%d\n",totsteps,ij,s,nbest,birth,newgene,nb1i,kch);
                }
                else {                                                      // use positional information to select ancestor (leave birth on)
                    nbest = s;
                    ancestor = newgene = parentid = 0ull;kch = 0;           // ancestor, newgene, parentid, kch initialized here to avoid warning below (not needed though)
                    for (nbmask=0ull,k=0;k<s;k++) nbmask |= 0x1ull<<((nb1i>>(k<<2))&0x7);    // check whether this needed here DEBUG !!!
                    if (!(nbest>1)) {                                       // for s==1 we define newgene ancestor immediately (s==0 does not reach here)
                        kch = nb1i&0x7;
                        ancestor = newgene = golg[nb[kch]];
                        parentid = golb[nb[kch]];
                    }
                    //if(ij==IJDEBUG) fprintf(stderr,"DEBUG Difference 2 at totsteps=%d ij=%d nbest=%d newgene=%llx nb1i=%llx kch=%d\n",totsteps,ij,nbest,newgene,nb1i,kch);
                }
                if(nbest>1 ) {
                    kch=selectdifft(nbest,nbmask,&crot,&kodd,&nsame);       // kch is chosen nb in range 0-7, nsame gives the number of undistinguished positions in canonical rotation
                    if(nsame) {
                        newgene = disambiguate(&kch, nb1i, nb, gol, golg, golb, nsame, &birth, &parentid, &ancestor, ij); // restore symmetry via one of 8 repscheme options
                        if(birth) statflag |= F_disambig;
                    }
                    else {
                        ancestor = newgene = golg[nb[kch]];
                        parentid = golb[nb[kch]];
                    }
                    //if(ij==IJDEBUG) fprintf(stderr,"DEBUG Difference 3 at totsteps=%d ij=%d nbest=%d newgene=%llx nb1i=%llx kch=%d\n",totsteps,ij,nbest,newgene,nb1i,kch);
                }
            }
            if (birth) {                                                    // ask again because disambiguate may turn off birth
                statflag |= F_birth;
                r2=1ull;                                                    // compute random events for single bit mutation, as well as mutation position nmut
                while (r2) {
                    RAND128P(randnr);                                       // inline exp so compiler recognizes auto-vec,
                    r2 = randprob(pmutmask,(unsigned int) randnr);
                    nmut = (randnr >> 56) & 0x3f;                           // choose mutation position for length 64 gene : from bits 56:61 of randnr
                    newgene = newgene ^ (r2<<nmut);                         // introduce single mutation with probability pmut = probmut
                    if(r2) {
                        statflag = statflag | F_mutation;
                        statflag = statflag | F_survmut;
                    }                }
                newgol[ij]  =  1ull;                                        // new game of life cell value: alive
                newgolg[ij] =  newgene;                                     // if birth then newgene
                newgolr[ij] = (golr[nb[kch]]<<4) | kch | 0x8ull;            // register ancestor offset index in record_of_dynamics_gene golr along with 0x8 for birth event
                if(diagnostics & diag_hash_genes) {
                    if(golij) hashdeletegene(golg[ij],golb[ij],"step %d hash delete error 1 in update_lut_sum, gene %llx not stored\n");
                    hashaddgene(ij,newgene,ancestor,newgolb+ij,parentid,statflag & F_mutation);
                }
                if (parentdies) {
                    unsigned int ij1 = nb[kch];
                    newgolgstats[ij1] = newgolgstats[ij1] | F_parent;
                }
                //if(ij==IJDEBUG) fprintf(stderr,"DEBUG Difference 4 at totsteps=%d ij=%d newgene=%llx nb1i=%llx kch=%d\n",totsteps,ij,newgene,nb1i,kch);
            }
        }
        if(!birth) {                                                       // need instead of else because if(birth) section may change value of birth
            if(golij) {                                                    // death/survival
                if(survive) {                                              // survival coded
                    statflag |= F_survival;
                    if (golgstats[ij]&F_survmut) statflag |= F_survmut;    // gene is non-replicated mutant survivor
                    newgol[ij]  = golij;                                   // new game of life cell value same as old
                    newgolg[ij] = golg[ij];                                // gene stays same
                    newgolb[ij] = golb[ij];
                    newgolr[ij] = (golr[ij]<<4) | (s-1);                   // register s-1 in record_of_dynamics_gene golr along with 0 value of bit3 for survival event
                }
                else {                                                     // death
                    statflag |= F_death;
                    newgol[ij]  = 0ull;                                    // new game of life cell value dead
                    newgolg[ij] = 0ull;                                    // gene dies
                    newgolb[ij] = 0ull;
                    newgolr[ij] = 0ull;
                    if(diagnostics & diag_hash_genes)
                        hashdeletegene(golg[ij],golb[ij],"step %d hash delete error 2 in update, gene %llx not stored\n");
                }
                //if(ij==IJDEBUG) fprintf(stderr,"DEBUG Difference 5 at totsteps=%d ij=%d nb1i=%llx\n",totsteps,ij,nb1i);
            }
            else {                                                         // empty and no birth, stays empty
                newgol[ij]  = golij;
                newgolg[ij] = golg[ij];
                newgolb[ij] = golb[ij];
                newgolr[ij] = golr[ij];
                //if(ij==IJDEBUG) fprintf(stderr,"DEBUG Difference 6 at totsteps=%d ij=%d nb1i=%llx\n",totsteps,ij,nb1i);
            }
        }
        if(newgol[ij]!=gols) {
            statflag |= F_notgolrul;
            if(newgol[ij]) statflag |= F_nongolchg;
        }
        if(golij) statflag |= F_golstate;                                   // this is the last gol state, not the new state
        if (parentdies) newgolgstats[ij] = newgolgstats[ij] | statflag;     // newgolgstats may already contain updated parenthood info F_parent
        else newgolgstats[ij] = statflag;
    
        //if(ij==IJDEBUG) fprintf(stderr,"DEBUG First error at totsteps=%d ij=%d s=%d newgolg=%llx newgolr=%llx\n",totsteps,ij,s,newgolg[ij],newgolr[ij]);
}
//........................................................................................................................................................
void extern inline finish_update(uint64_t newgol[], uint64_t newgolg[],uint64_t newgolgstats[],uint64_t newgolb[], uint64_t newgolr[], int nbshist[]) {
    int ij,k;
    uint64_t statflag;
    if(parentdies) {
        for (ij=0; ij<N2; ij++) {
            statflag = newgolgstats[ij];
            if(gol[ij] && !(statflag&F_birth) && !(statflag&F_death) && (statflag&F_parent)) {
                newgol[ij]  = 0ull;                                    // new game of life cell value dead
                newgolg[ij] = 0ull;                                    // gene dies
                newgolb[ij] = 0ull;
                newgolr[ij] = 0ull;
                if(diagnostics & diag_hash_genes)
                    hashdeletegene(golg[ij],golb[ij],"step %d hash delete error 2 in update, gene %llx not stored\n");
                newgolgstats[ij] |= F_parentaldeath;
            }
        }
    }
    for (ambigsum=0,ij=0; ij<N2; ij++) ambigsum += (newgolgstats[ij]&F_disambig) ? 1 : 0;
    
    for(k=0;k<8;k++) nbshist[k]=0;
    for (ij=0; ij<N2; ij++) {
        if(newgolgstats[ij]&F_birth)
            for (k=0; k<8; k++) nbshist[k] += (newgolgstats[ij]&(0x1ull<<(16+k))) ? 1: 0;
    }
    
    if(randominflux) random_influx(newgol,newgolg,newgolb,newgolr);
    if(vscrolling) v_scroll(newgol,newgolg,newgolb,newgolr);
    if ((colorfunction == 8) || (colorfunction2 == 8)) packandcompare(newgol,working,golmix);
    if(diagnostics & diag_component_labels) ncomponents=extract_components(newgol);
    if(diagnostics & diag_hash_genes) {
        for (ij=0; ij<N2; ij++) {       // complete missing hash table records of extinction and activities
            if(gol[ij]) hashgeneextinction(golg[ij],"hash extinction storage error %d in update at step %d, gene %llx not stored\n");     // [**gol**]
            if(newgol[ij]) hashgeneactivity(newgolg[ij],"hash activity storage error in update, gene %llx not stored\n");   // [**gol**]
        }
    }
    // if(diagnostics & diag_hash_patterns) qimage = quadimage(newgol,&patt,log2N); // quadtree hash of entire image
}
//---------------------------------------------------------------- update_lut_sum -----------------------------------------------------------------------
void update_lut_sum(uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[],uint64_t golr[],uint64_t newgol[], uint64_t newgolg[], uint64_t newgolgstats[], uint64_t newgolb[],uint64_t newgolr[]){    // selection models 8,9
// this version should work even if extra information is packed in the higher bits of gol: previous version relabelled to update_lut_sumx
// update GoL for toroidal field which has side length which is a binary power of 2
// encode without if structures for optimal vector treatment
/*
        0 <= s <= 8
        genome =
        8 bits for each possible s value for birth (center site 0) : exclude s=0 no birth of isolated live cells
        8 bits for each possible s value for survival/death (center site 1) : exclude s=0 no survival of isolated cells
        1. fixed length encoding for selection=8: using b0 for s=1, b1 for s=2, etc
        s-1   76543210    76543210
        GOL = 00000100(0)|00000110(0)
            = 0000 0100 0000 0110
            = 0x0406
        2. variable length encoding for selection=9: using 4 bit patterns anywhere on 4-bit raster in genome
        note that one entry suffices for any lut rule: duplicate rules may be used to encode robustness
        of four bits 0-3, bit 3 encodes birth/survival as 1/0. Bits 0-2 encode the value of s (0-7, s=8 excluded).
        GOL = 0xb32 or 0x23b2b for example (up to 64 bits)
                                                                                                                */
    int s, s1, s2or3, k, kch;
    unsigned int rulemodij;
    int nb[8],  ij, i, j, jp1, jm1, ip1, im1;
    uint64_t genecode, genecodeb, genecodes, gols, golij, nb1i, nbmask, found;
    uint64_t  survive, birth, overwrite, survivalgene, smask, bmask, ncodingmask, allcoding;

    canonical = repscheme & R_2_canonical_nb;                                   // set global choice of canonical rotation bit choice for selectdifftx
    survivalgene = repscheme & R_0_survivalgene;                                // gene determining survival is 1: central gene 0: determined by neighbours
    smask = (uint64_t) survivalmask;                                            // convert to 64 bit mask for efficient usage here
    bmask = (uint64_t) birthmask;
    parentdies = (repscheme & R_3_parentdies) ? 1 : 0;
    if(parentdies) for (ij=0; ij<N2; ij++) newgolgstats[ij] = 0ull;             // need to update statistics of neighbours with parenting information, so init required
    
    ncodingmask = (1ull<<ncoding)-1ull;                                         // mask for number of bits coding for each lut rule: <=4 for birth and survival case
    if (ncoding==4)         allcoding = 0xffffffffffffffff;                     // mask for total gene coding region
    else if (ncoding ==2)   allcoding = 0xffffffff;
    else                    allcoding = 0xffff;
    
    /* DEBUG if(totsteps<5 && selection==9) {
        retrievegols(totsteps,newgol,newgolg,newgolgstats,newgolb,newgolr);     // previous simulation data placed in newgol...
        for (ij=0; ij<N2; ij++) {
            if(gol[ij]!=newgol[ij]) {
                fprintf(stderr,"NEW step %d gol difference at ij = %d\n",totsteps,ij);
            }
        }
    }*/
  
    for (ij=0; ij<N2; ij++) {                                                   // loop over all sites of 2D torus with side length N
        i = ij & Nmask;  j = ij >> log2N;                                       // row & column
        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                       // toroidal (j+1)*N and (j-1)*N
        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                             // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;               // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (s=0,nb1i=0ull,nbmask=0ull,k=0;k<8;k++) {                           // packs non-zero nb indices in first up to 8*4 bits
            gols=gol[nb[k]];                                                    // whether neighbor is alive
            s += gols;                                                          // s is number of live nbs
            nb1i = (nb1i << (gols<<2)) + (gols*k);                              // nb1i is packed list of live neighbour indices (each in range 0-7)
            nbmask |= (gols << k);
        }
        survive = birth = 0ull;
        golij=gol[ij];
        if (s) {
            s1=s-1;
            s2or3 = (s>>2) ? 0ull : (s>>1);                                     // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
            gols = s2or3 ? (golij ? 1ull : (s&1ull ? 1ull : 0ull )) : 0ull;     // GoL calculation next state for non-genetic gol plane
            rulemodij = (rulemod&0x2) ? (ij>=(N2>>1) ? 1 : 0) : (rulemod&0x1);  // if rulemod bit 1 is on then split into half planes with/without mod
            rulemodij = (rulemod&0x4) ? membrane : (rulemod&0x1);               // if rulemod bit 2 then activate membrane of death overruling previous line
            if(rulemodij==1) {
                overwrite = overwritemask&(0x1ull<<s1);
                if (selection==9) {                                             // selection == 9 : NB selection 9 does not use ncoding to determine coding length
                    for (genecodes=genecodeb=(repscheme&R_1_nb_OR_AND)?0ull:0xffff,k=0;k<s;k++) {// decodes genes using variable position encoding: fixed length genecode from nbs variable position encoding
                        kch = (nb1i>>(k<<2))&0x7;
                        if (!survivalgene && golij) {
                            PATTERN4(golg[nb[kch]], (s&0x7), found);          //survival?
                            if(repscheme&R_1_nb_OR_AND)
                                genecodes |= found? (1ull << s1) : 0ull;
                            else
                                genecodes &= found? (1ull << s1) : 0ull;
                        }
                        if (overwrite || !golij) {
                            PATTERN4(golg[nb[kch]], (s|0x8), found);          //birth?
                            if(repscheme&R_1_nb_OR_AND)
                                genecodeb |= found ? (1ull << (s1+8)) : 0ull;
                            else {
                                genecodeb &= found ? (1ull << (s1+8)) : 0ull;
                                //if(ij==106195) fprintf(stderr,"DEBUG Genecode sel 9 at totsteps=%d ij=%d k=%d golg=%llx found=%llx genecode=%llx survivalgene=%llx\n",totsteps,ij,k,golg[nb[k]],found,genecodeb,survivalgene);
                            }
                        }
                        // if(ij==106195) fprintf(stderr,"DEBUG Genecode sel 9 at totsteps=%d ij=%d k=%d golg=%llx golij=%llx\n",totsteps,ij,k,golg[nb[k]],golij);
                    }
                    if (survivalgene && golij) {                                // survival determined by central gene in this case
                                PATTERN4(golg[ij], (s&0x7), found);             // survival?
                                genecodes =  found? 1ull << s1 : 0ull;
                    }
                    genecode=genecodes|genecodeb;
                    //if(ij==106195) fprintf(stderr,"DEBUG Genecode sel 9 at totsteps=%d ij=%d genecode=%llx\n",totsteps,ij,genecode);
                    survive = ((genecode&smask)>>s1) & 0x1ull;
                    genecode>>=8;
                    if (overwrite || !golij) birth   = ((genecode&bmask)>>s1) & 0x1ull;
                }
                else {                                                          // selection == 8
                    if(repscheme&R_1_nb_OR_AND)
                        for (genecode=0ull,k=0;k<8;k++) {                       // decodes genes with fixed length encoding by OR
                            kch = (nb1i>>(k<<2))&0x7;
                            // genecode |= ((gol[nb[kch]]&0x1ull)?golg[nb[kch]]:0ull);   // OR of live neighbours encodes birth rule & survival rule
                            genecode |= golg[nb[kch]];   // OR of live neighbours encodes birth rule & survival rule
                        }
                    else
                        for (genecode=allcoding,k=0;k<s;k++) {                   // decodes genes with fixed length encoding by AND
                            kch = (nb1i>>(k<<2))&0x7;
                            // genecode &= ((gol[nb[kch]]&0x1ull)?golg[nb[kch]]:allcoding); // AND of live neighbours encodes birth rule & survival rule
                            genecode &= golg[nb[kch]];                           // AND of live neighbours encodes birth rule & survival rule
                        }
                    if(survivalgene) genecode = (genecode&(0xffffffffull<<32)) | (golg[ij]&(8ull*ncoding-1ull));   // if central gene codetermines survival DEBUG CHECK
                    //if(ij==106195) fprintf(stderr,"DEBUG Genecode sel 8 at totsteps=%d ij=%d genecode=%llx\n",totsteps,ij,genecode);
                    //if(ij==106195) fprintf(stderr,"DEBUG Genecode sel 8 at totsteps=%d ij=%d genecode=%llx golg1=%llx golg2=%llx\n",totsteps,ij,genecode,golg[nb[nb1i&0x7]],golg[nb[(nb1i>>(1<<2))&0x7]]);
                    survive=(((genecode>>((s-1)*ncoding)) & ncodingmask) == ncodingmask) && ((smask>>s1)&1ull) ? 1ull : 0ull;
                    if (overwrite || !golij) birth=(((genecode>>((8+s1)*ncoding)) & ncodingmask) == ncodingmask) && ((bmask>>s1)&1ull) ? 1ull : 0ull;
                }
            }
            else if (rulemodij==2){                                          // hard death on membrane defined above via macro "membrane"
                survive = 0ull;
                birth = 0ull;
            }
            else {
                    survive = s2or3;
                    birth = s2or3&s&0x1ull&~golij;
            }
        }
        else survive = birth = s2or3 = gols = 0ull;
        finish_update_ij(ij,s,golij,gols,nb1i,nbmask,nb,survive,birth,gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);
    }  // end for ij

    finish_update(newgol, newgolg, newgolgstats, newgolb, newgolr, nbshist);
    
    /* if(totsteps<5 && selection == 8) {
        savegols(totsteps,gol,golg,golgstats,golb,golr);
    } */

}
//---------------------------------------------------------------- update_lut_dist ----------------------------------------------------------------------
void update_lut_dist(uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[],uint64_t golr[],uint64_t newgol[], uint64_t newgolg[], uint64_t newgolgstats[], uint64_t newgolb[],uint64_t newgolr[]) {     // selection models 10,11
// update GoL for toroidal field which has side length which is a binary power of 2
// encode without if structures for optimal vector treatment
/*
    configurations are distinguished by number of ones in two classes (corner,edge-centred) of the 8-bit live neighbour pattern
    for the different s values                 0  1  2  3  4  5  6  7  8
    there are a nr n of partitions             1  2  3  4  5  4  3  2  1   total 25   with se the number in edge centred sites NSEW
    chosen nr of partitions for exploration    0  2  3  4  5  4  3  2  0   total 23
    if we exclude the cases s = 0,8, then there are 23 bits required to distinguish these cases
    1. Fixed length encdoing has 1 bit per LUT entry: 23 survival bits and then 23 birth bits
    2. Varibale length encoding uses two bytes : upper byte for s-1 and lower byte for bitmask of which se 0-3 are on
                                                                                                                */
    int s, s1, se, s0, s2or3, k, kch;
    unsigned int rulemodij;
    int nb[8], ij, i, j, jp1, jm1, ip1, im1;
    uint64_t genecode, survivenbs, birthnbs, gols, golij, nb1i, nbmask, found;
    uint64_t survive, birth, overwrite, survivalgene, smask, bmask, gene;
    
    static uint64_t summasks[7] = {0x3ull,0x7ull,0xfull,0x1full,0xfull,0x7ull,0x3ull};
    static int sumoffs[7] = {0,2,5,9,14,18,21};
    // static int first = 1;

    canonical = repscheme & R_2_canonical_nb;                                      // set global choice of canonical rotation bit choice for selectdifftx
    survivalgene = repscheme & R_0_survivalgene;                                   // gene determining survival is 1: central gene 2: determined by neighbours
    smask = (uint64_t) survivalmask;                                               // convert to 64 bit mask for efficient usage here
    bmask = (uint64_t) birthmask;
    parentdies = (repscheme & R_3_parentdies) ? 1 : 0;
    if(parentdies) for (ij=0; ij<N2; ij++) newgolgstats[ij] = 0ull;                // need to update statistics of neighbours with parenting information, so init required

    /* if(totsteps<5) {
        retrievegols(totsteps,newgol,newgolg,newgolgstats,newgolb,newgolr);     // previous simulation data placed in newgol...
        for (ij=0; ij<N2; ij++) {
            if(gol[ij]!=newgol[ij]) {
                if(first) {fprintf(stderr,"DEBUG NEW step %d gol difference at ij = %d this gol %llx file gol %llx\n",totsteps,ij,gol[ij],newgol[ij]);first=0;}
            }
        }
    } */
    
    for (ij=0; ij<N2; ij++) {                                                      // loop over all sites of 2D torus with side length N
        i = ij & Nmask;  j = ij >> log2N;                                          // row & column
        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                          // toroidal (j+1)*N and (j-1)*N
        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                                // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;                  // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (s=se=0,nb1i=0ull,nbmask=0ull,k=0;k<8;k++) {                           // packs non-zero nb indices in first up to 8*4 bits
            gols=gol[nb[k]];                                                       // whether neighbor is alive
            s += gols;                                                             // s is number of live nbs
            se += k&0x1&gols;                                                      // se is number of edge-centred live neighbours (odd k)
            nb1i = (nb1i << (gols<<2)) + (gols*k);                                 // nb1i is packed list of live neighbour indices
            nbmask |= (gols << k);
        }
        survive = birth = 0ull;
        golij = gol[ij];
        
        if (s>0 && s<8) {
            s1 = s-1;
            s2or3 = (s>>2) ? 0ull : (s>>1);                                        // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
            gols = s2or3 ? (golij ? 1ull : (s&1ull ? 1ull : 0ull )) : 0ull;        // GoL calculation next state for non-genetic gol plane
            rulemodij = (rulemod&0x2) ? (ij>=(N2>>1) ? 1 : 0) : (rulemod&0x1);     // if rulemod bit 1 is on then split into half planes with/without mod
            rulemodij = (rulemod&0x4) ? membrane : (rulemod&0x1);                  // if rulemod bit 2 then activate membrane of death
            if (rulemodij==1) {                // NB need to put gene calculation outside so that we can do genetic propagation with GoL rulemod off
                overwrite = overwritemask&(0x1ull<<s1);
                overwrite = (overwrite || !golij) ? 1ull : 0ull;
                if (golij) survive = (smask>>sumoffs[s1])&summasks[s1] ? 1ull : 0ull;
                if (overwrite) birth = (bmask>>sumoffs[s1])&summasks[s1] ? 1ull : 0ull;
                if (survive|birth) {
                    s0 = (s-4 > 0 ? s-4 : 0);
                    if (golij)   survive = (smask>>(sumoffs[s1]+se-s0))&0x1ull;     // refine decisions for specific combination of s and se
                    if (overwrite) birth = (bmask>>(sumoffs[s1]+se-s0))&0x1ull;     // only allowed if birth,survivalmask permits (ask this before consulting genes)
                    if (survive|birth) {                                            // complete determination of birth or survival
                        if (selection==11) {
                            survivenbs=birthnbs=(repscheme&R_1_nb_OR_AND) ? 0ull : 1ull; // initialization for OR or AND of neighbors resp.
                            for (k=0;k<s;k++)  {                                    // decodes genes with variable position encoding only for current s,se
                                kch = (nb1i>>(k<<2))&0x7;
                                gene = golg[nb[kch]] & (0xf0f0f0f0f0f0f0f0ull | (0x0101010101010101ull<<(se-s0))); // focus gene down to specific required se bit for exact matching
                                if (!survivalgene && golij) {                 // coding is 8 bits [(b/s) (s1 2 1 0) (se subset mask 3 2 1 0)]  (exception case s==4,se==4 is x0000001)
                                    PATTERN8(gene, (s==4 && se==4) ? 0x01 : ((s<<4)|(1<<(se-s0))), found);    //survival rule found? final decision for survival
                                    if(repscheme&R_1_nb_OR_AND)
                                        survivenbs |= found? 1ull : 0ull;           // incl. special case codes for 5th case se==4 for s==4
                                    else
                                        survivenbs &= found? 1ull : 0ull;           // incl. special case codes for 5th case se==4 for s==4
                                }
                                if (overwrite) {
                                    PATTERN8(gene, (s==4 && se==4) ? 0x81 :(((8|s)<<4)|(1<<(se-s0))), found); //birth rule found? final decision for birth
                                    if(repscheme&R_1_nb_OR_AND)
                                        birthnbs |= found? 1ull : 0ull;
                                    else
                                        birthnbs &= found? 1ull : 0ull;
                                }
                                //if(ij==IJDEBUG) fprintf(stderr,"DEBUG Genecode sel 11 at totsteps=%d ij=%d k=%d golg=%llx golij=%llx\n",totsteps,ij,k,golg[nb[k]],golij);
                            }
                            if (survivalgene && golij) {                            // survival determined by central gene in this case
                                    PATTERN8(golg[ij], (s==4 && se==4) ? 0x01 : ((s<<4)|(1ull<<(se-s0))), found);   // survival rule found? final decision for survival
                                    survivenbs |= found? 1ull : 0ull;               // special case codes for 5th case se==4 for s==4
                            }
                            birth &= birthnbs;
                            survive &= survivenbs;
                        }
                        else {                                                      // selection == 10
                            if(repscheme & R_1_nb_OR_AND)
                                for (genecode=0ull,k=0;k<s;k++) {                   // decodes genes with fixed length encoding by OR
                                    kch = (nb1i>>(k<<2))&0x7;
                                    genecode |= golg[nb[kch]];                      // OR of live neighbours encodes birth rule & survival rule
                                }
                            else
                                for (genecode=~0ull,k=0;k<s;k++) {                  // decodes genes with fixed length encoding by AND
                                    kch = (nb1i>>(k<<2))&0x7;
                                    genecode &= golg[nb[kch]];                      // AND of live neighbours encodes birth rule & survival rule
                                }
                            if(survivalgene) genecode = (genecode&(0xffffffffull<32)) | (golg[ij]&0xffffffffull);     // if central gene determines survival
                            genecode &= (bmask << 32)|smask;                        // just to be clear, not required as now already tested above
                            //if(ij==IJDEBUG) fprintf(stderr,"DEBUG Genecode sel 10 at totsteps=%d ij=%d genecode=%llx\n",totsteps,ij,genecode);
                            //if(ij==IJDEBUG) fprintf(stderr,"DEBUG Genecode sel 10 at totsteps=%d ij=%d genecode=%llx golg1=%llx golg2=%llx\n",totsteps,ij,genecode,golg[nb[nb1i&0x7]],golg[nb[(nb1i>>(1<<2))&0x7]]);
                            if (golij) survive &= (genecode>>(sumoffs[s1]+(se-s0)))&0x1ull;        // 28.7.2019 changed from |=
                            if (overwrite) birth &= (genecode>>(32+sumoffs[s1]+(se-s0)))&0x1ull;
                        }
                    }
                }
            }
            else if (rulemodij==2){                                                 // hard death on membrane defined above via macro "membrane"
                survive = 0ull;
                birth = 0ull;
            }
            else {                                                                  // GOL rule
                survive = s2or3;
                birth = s2or3&s&0x1&~golij;
            }
            // if(s>2)  survive=0ull;
        }
        else survive = birth = s2or3 = gols = 0ull;

        finish_update_ij(ij,s,golij,gols,nb1i,nbmask,nb,survive,birth,gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);
    }  // end for ij

    finish_update(newgol, newgolg, newgolgstats, newgolb, newgolr, nbshist);
}
//---------------------------------------------------------------- update_lut_canon_rot -----------------------------------------------------------------
void update_lut_canon_rot(uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[],uint64_t golr[],uint64_t newgol[], uint64_t newgolg[], uint64_t newgolgstats[], uint64_t newgolb[],uint64_t newgolr[]) {     // selection models 12,13
/*
    configurations are distinguished by number of ones and the canonical rotation of the 8-bit live neighbour pattern
    the canonical rotation is the rotation with the minimum numerical value as 8-bit number (bits are numbered clockwise from 0-7)
    for the different s values                 0  1  2  3  4  5  6  7  8
    there are a nr of canonical rotations      1  1  4  7 10  7  4  1  1   total 36
    if we exclude the cases s=0,1,7,8 which can create growth artefacts, then there are 32 bits
 
    1. Fixed length encoding: Survival: s=2 bits 0-3, s=3 bits 4-10, s=4 bits 11-20, s=5 bits 21-27, s=6 bits 28-31; Birth: same + 32
    2. Modular encoding: Up to 6 10-bit modules. Lower 8 bits stored in first 6*8=48 bits. Upper 2-bits stored pairwise in 12 bits 48-59.
       Module: 1     0    >=48 ...  7     6     5     4     3     2     1     0       we abbreviate possible values of crot mod 5 as cr0-4 in this description
               cr4   cr3            s3    s2    s1    s0    crot5 cr2   cr1   cr0     8-bit pattern matching via remapping cr0-4 to lowest 3 bits: 0->001 1->010 2->100 3->001 4->010
                                    s3    s2    s1    s0    crot5 0/cr2 cr4/1 cr3/0   remapping of bits depending on whether crot mod 5 is >2 / <=2
                                                                                                */
    int s, smid, s2, s2or3, k, k1, kodd, crot, crot5, crotmod5, pat, found, nsame;
    uint64_t survive,birth,survivenbs,birthnbs,overwrite,survivalgene,smask,bmask,rulemodij,golij;
    static uint64_t summasks[5] = {0xfull,0x7full,0x3ffull,0x7full,0xfull};
    static int sumoffs[5] = {0,4,11,21,28};
    int nb[8], ij, i, j, jp1, jm1, ip1, im1;
    unsigned int kch=0;
    uint64_t genecode, genecode1, gols, nb1i, nbmask;
    uint64_t newgene;
    uint64_t statflag;

    canonical = repscheme & R_2_canonical_nb;                                       // set global choice of canonical rotation bit choice for selectdifftx
    survivalgene = repscheme & R_0_survivalgene ? 1ull : 0ull;                      // gene determining survival is 1: central gene 2: determined by neighbours
    smask = (uint64_t) survivalmask;                                                // 32 bits of survivalmask used to limit space of rules, convert to 64 bit masks for efficient usage here
    bmask = (uint64_t) birthmask;                                                   // 32 bits of birthmask
    parentdies = (repscheme & R_3_parentdies) ? 1 : 0;
    
    if(parentdies) for (ij=0; ij<N2; ij++) newgolgstats[ij] = 0ull;                 // need to update statistics of neighbours with parenting information, so init required

    for (ij=0; ij<N2; ij++) {                                                       // loop over all sites of 2D torus with side length N
        i = ij & Nmask;  j = ij >> log2N;                                           // row & column
        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                           // toroidal (j+1)*N and (j-1)*N
        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                                 // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;                   // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (s=0,nb1i=0ull,nbmask=0ull,k=0;k<8;k++) {                               // packs non-zero nb indices in first up to 8*4 bits
            gols=gol[nb[k]];                                                        // whether neighbor is alive
            s += gols;                                                              // s is number of live nbs
            nb1i = (nb1i << (gols<<2)) + (gols*k);                                  // nb1i is packed list of live neighbour indices
            nbmask |= (gols << k);
        }
        smid = s>1 && s<7; s2 = s-2;                                                // s in mid-range for possible lut rule
        statflag = F_livenbs & (nbmask<<16);                                        // record live neighbour pattern
        statflag |= F_s_live & (s&0x7);                                             // requires F_s_live to be in lowest 3 bits : val 8 mapped to 0
        newgene = survive = birth = 0ull;
        golij = gol[ij];

        if (smid) {
            s2or3 = (s>>2) ? 0ull : (s>>1);                                         // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
            gols = s2or3 ? (gol[ij] ? 1ull : (s&1ull ? 1ull : 0ull )) : 0ull;       // GoL calculation next state for non-genetic gol plane
            rulemodij = (rulemod&0x2) ? (ij>=(N2>>1) ? 1 : 0) : (rulemod&0x1);      // if rulemod bit 1 is on then split into half planes with/without mod
            rulemodij = (rulemod&0x4) ? membrane : (rulemod&0x1);                   // if rulemod bit 2 then activate membrane of death
            if(rulemodij==1) {
                overwrite = s ? (overwritemask>>(s-1))&0x1ull : 0ull;               // allow birth to overwrite occupied cell = survival in GoL
                overwrite = overwrite | (~gol[ij] & 0x1ull);                        // either central cell is empty or overwrite bit set is required for birth
                if (gol[ij]) survive = (smask>>sumoffs[s2])&summasks[s2] ? 1ull : 0ull;
                if (overwrite) birth = (bmask>>sumoffs[s2])&summasks[s2] ? 1ull : 0ull;
                if (survive|birth) {
                    for (k=0,nbmask=0;k<8;k++) nbmask |= (gol[nb[k]]<<k);           // constuct mask of live bits for 8 neighbours
                    kch = selectdifft(s, nbmask, &crot, &kodd, &nsame);             // find the canonical rotation index of the live neighbour configuration
                    survive&= (smask>>(sumoffs[s2]+crot))&0x1ull;                   // refine decisions for specific canonical rotation configuration
                    birth  &= (bmask>>(sumoffs[s2]+crot))&0x1ull;                   // only allowed if birth/survivemask permits (ask this before consulting genes)
                    if (survive|birth) {                                            // complete determination of birth or survival
                        if (selection==13) {                                        // modular gene encoding one of up to two 5-bit subsets of crot (e.g. for s=3,4 5+2 and 5+5 resp.)
                            survivenbs=birthnbs=(repscheme&R_1_nb_OR_AND) ? 0ull : 1ull; // initialization for OR or AND of neighbors resp.
                            crot5 = crot > 4 ? 1 : 0;                               // crot is in range 0-9
                            crotmod5 = crot-5*crot5;                                // the 5 bits indexed by crotmod5 are stored in two places 0-2 in 8-bit word and 3-4 in pairs bits 48+
                            pat = (s<<4)|(crot5<<3)|(0x1<<((crotmod5>2)?crotmod5-3:crotmod5));
                            genecode = 0ull;
                            for (k=0;k<8;k++) {                                     // decodes genes with variable position encoding only for current s,crot
                                if (gol[nb[k]]) {                                   // combine information from genes of all live neighbours
                                    if ((~survivalgene & gol[ij]) | overwrite) {
                                        genecode = golg[nb[k]];
                                        if(crotmod5 > 2) for (genecode1=0ull,k1=0;k1<6;k1++) genecode1 |= ( (genecode>>(48+(k1<<1))) & (0x1ull<<(crotmod5-3)) ) << (k1<<3);
                                        else genecode1 = genecode & (0x010101010101 << crotmod5); // prepare lookup in all 6 modules on gene : in this case crot subset bits in right place
                                        genecode &= 0xf8f8f8f8f8f8;
                                        genecode |= genecode1;                      // replace lowest 3 bits of 8-bit part of 10-bit modules with appropriate crot subset bits
                                    }
                                    if (~survivalgene & gol[ij]) {
                                        PATTERN8(genecode, pat, found);
                                        if(repscheme&R_1_nb_OR_AND)
                                            survivenbs |= found? 1ull : 0ull;
                                        else
                                            survivenbs &= found? 1ull : 0ull;
                                    }
                                    if (overwrite) {
                                        PATTERN8(genecode, (0x80|pat), found);
                                        if(repscheme&R_1_nb_OR_AND)
                                            birthnbs |= found? 1ull : 0ull;
                                        else
                                            birthnbs &= found? 1ull : 0ull;
                                    }
                                }
                            }
                            /* fprintf(stderr,"lut_canon step %d ij %d s %d gol[ij] %llx crot %d (%d,%d) pat %x genecode %llx overwrite %llx survive %llx birth %llx\n",
                                                     totsteps,ij,s,gol[ij],crot,crot5,crotmod5,pat,genecode,overwrite,survive,birth); */
                            // if (overwrite && !gol[ij] && s==3) fprintf(stderr,"overwrite and golij!!!!!!!!!!!!!\n");
                            if (survivalgene & gol[ij]) {                           // survival determined by central gene in this case
                                genecode = golg[ij];
                                if(crotmod5 > 2) for (genecode1=0ull,k1=0;k1<6;k1++) genecode1 |= ((genecode>>(48+(k1<<1)))&(0x1<<(crotmod5-3)))<<(k1<<3);
                                else genecode1 = genecode & (0x010101010101 << crotmod5);
                                genecode &= 0xf8f8f8f8f8f8;
                                genecode |= genecode1;
                                PATTERN8(genecode, pat, found);
                                survivenbs = found? 1ull : 0ull;
                            }
                            birth &= birthnbs;
                            survive &= survivenbs;
                        }
                        else {                                                      // selection == 12
                            if(repscheme&R_1_nb_OR_AND)
                                for (genecode=0ull,k=0;k<8;k++)                     // decodes genes with fixed length encoding by OR
                                    genecode |= (gol[nb[k]]?golg[nb[k]]:0ull);      // OR of live neighbours encodes birth rule & survival rule
                            else
                                for (genecode=~0ull,k=0;k<8;k++)                    // decodes genes with fixed length encoding by AND
                                    genecode &= (gol[nb[k]]?golg[nb[k]]:~0ull);     // AND of live neighbours encodes birth rule & survival rule
                            if(survivalgene) genecode = (genecode&(0xffffffffull<<32)) | (golg[ij]&0xffffffffull);     // if central gene determines survival
                            genecode&=(bmask<<32)|smask;                            // actually no longer needed since test done above
                            if (gol[ij]) survive &= (genecode>>(sumoffs[s2]+crot))&0x1ull;         // 28.7.2019 changed from |=
                            if (overwrite) birth &= (genecode>>(32+sumoffs[s2]+crot))&0x1ull;
                        }
                    }
                }
            }
            else if (rulemodij==2){                                                 // hard death on membrane defined above via macro "membrane"
                survive = 0ull;
                birth = 0ull;
            }
            else {
                survive = s2or3;
                birth = s2or3&s&0x1&~gol[ij];
            }
        }
        else survive = birth = s2or3 = gols = 0ull;
        
        finish_update_ij(ij,s,golij,gols,nb1i,nbmask,nb,survive,birth,gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);
    }  // end for ij

    finish_update(newgol, newgolg, newgolgstats, newgolb, newgolr, nbshist);
}
//---------------------------------------------------------------- update_lut_2Dsym ---------------------------------------------------------------------
void update_lut_2D_sym(uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[],uint64_t golr[],uint64_t newgol[], uint64_t newgolg[], uint64_t newgolgstats[], uint64_t newgolb[],uint64_t newgolr[]) {     // selection models 14,15
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
    1. Fixed length encoding: Survival: s=0 bit 0, s=1 bits 1,2, s=2 bits 3-8, s=3 bits 9-18, s=4 bits 19-31; Birth: same + 32
    2. Modular encoding: Up to 5 12-bit modules. Lower 8 bits in up to 5 modules in bits 0-39, upper 4 bits in up to 5 modules in bits 40-59

                                                                                                                */
    int s, se, slow, s2or3, k, k1, kodd, coff, crot, found, coffdiv6, coffmod6, pat, nsame;
    uint64_t survive,birth,survivenbs,birthnbs,overwrite,survivalgene,smask,bmask,rulemodij,golij;
    static uint64_t summasks[5] = {0x1ull,0x3ull,0x3full,0x3ffull,0x1fffull};  // masks for s= 0,1,2,3,4 with nr cases 1,2,6,10,13
    static int sumoffs[9] = {0,1,3,9,19,32,42,48,50};                          // cumulative offsets to start of coding region for s = 0,1,2,3,4,5,6,7,8 - only s<5 used
    static int csumoffs[9] = {0,2,4,12,26,46,60,68,2};                         // start of indexing in confoffs for crot,kodd lookup for s = 0,1,2,3,4,5,6,7,8
    static unsigned char confoffs[2+2+8+14+20+14+8+2+2] = {0,0, 0,0, 0,0,1,2,3,3,4,5, 0,1,2,3,3,2,4,5,6,7,4,5,8,9, 0,0,1,2,3,4,1,2,5,6,7,8,9,9,10,10,8,7,11,12,
                                        0,1,2,3,3,2,4,5,6,7,4,5,8,9, 0,0,1,2,3,3,4,5, 0,0, 0,0}; // look up for crot*2+kodd to gene bit offset
    static unsigned char lut[256];
    int nb[8], ij, i, j, jp1, jm1, ip1, im1;
    unsigned int kch=0;
    uint64_t genecode, genecode1, gols, nb1i, nbmask;
    uint64_t newgene;
    uint64_t statflag;

    static int first = 1;

    if (first) {                                                                // call selectdifft once at start for all possible nbmask vals to accelerate processing
        first = 0;
        for (nbmask=0;nbmask<256;nbmask++) {
            POPCOUNT64C(nbmask,s);
            kch = selectdifft(s, nbmask, &crot, &kodd, &nsame);
            coff = confoffs[csumoffs[s]+(crot<<1)+kodd];
            lut[nbmask]=coff;
            // fprintf(stderr,"LUT nbmask %2llx s %1d goff %2d soff %2d coff %2d crot %2d kodd %1d nsame %2d\n",nbmask,s,sumoffs[s]+coff,sumoffs[s],coff,crot,kodd,nsame);
        }
    }

    canonical = repscheme & R_2_canonical_nb;                                  // set global choice of canonical rotation bit choice for selectdifftx
    survivalgene = repscheme & R_0_survivalgene;                               // gene determining survival is 1: central gene 2: determined by neighbours
    smask = (uint64_t) survivalmask;                                           // 32 bits of survivalmask used to limit space of rules, convert to 64 bit masks for efficient usage here
    bmask = (uint64_t) birthmask;                                              // 32 bits of birthmask
    parentdies = (repscheme & R_3_parentdies) ? 1 : 0;

    if(parentdies) for (ij=0; ij<N2; ij++) newgolgstats[ij] = 0ull;            // need to update statistics of neighbours with parenting information, so init required

    for (ij=0; ij<N2; ij++) {                                                  // loop over all sites of 2D torus with side length N
        i = ij & Nmask;  j = ij >> log2N;                                      // row & column
        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                      // toroidal (j+1)*N and (j-1)*N
        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                            // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;              // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (s=se=0,nb1i=0ull,nbmask=0ull,k=0;k<8;k++) {                       // packs non-zero nb indices in first up to 8*4 bits
            gols=gol[nb[k]];                                                   // whether neighbor is alive
            s += gols;                                                         // s is number of live nbs
            se += k&0x1&gols;                                                  // se is number of edge-centred live neighbours (odd k)
            nb1i = (nb1i << (gols<<2)) + (gols*k);                             // nb1i is packed list of live neighbour indices
            nbmask |= (gols << k);
        }
        slow = s<5;                                                            // s in low-range 0-4 for possible lut rule
        statflag = F_livenbs & (nbmask<<16);                                   // record live neighbour pattern
        statflag |= F_s_live & (s&0x7);                                        // requires F_s_live to be in lowest 3 bits : val 8 mapped to 0
        newgene = survive = birth = 0ull;
        golij = gol[ij];

        if (slow) {
            s2or3 = (s>>2) ? 0ull : (s>>1);                                    // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
            gols = s2or3 ? (gol[ij] ? 1ull : (s&1ull ? 1ull : 0ull )) : 0ull;  // GoL calculation next state for non-genetic gol plane
            rulemodij = (rulemod&0x2) ? (ij>=(N2>>1) ? 1 : 0) : (rulemod&0x1); // if rulemod bit 1 is on then split into half planes with/without mod
            rulemodij = (rulemod&0x4) ? membrane : (rulemod&0x1);              // if rulemod bit 2 then activate membrane of death
            if(rulemodij==1) {
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
                    survive= (smask>>(sumoffs[s]+coff))&0x1ull;                // refine decisions for specific canonical rotation configuration
                    birth  = (bmask>>(sumoffs[s]+coff))&0x1ull;                // only allowed if birth/survivemask permits (ask this before consulting genes)
                    if (survive|birth) {                                       // complete the determination of birth or survival
                        if (selection==15) {                                   // selection == 15 variable length encoding
                            survivenbs=birthnbs=(repscheme&R_1_nb_OR_AND) ? 0ull : 1ull; // initialization for OR or AND of neighbors resp.
                            for (k=0;k<8;k++) {                                // decodes genes with variable position encoding only for current s,crot
                                if (gol[nb[k]]) {                              // combine information from genes of all live neighbours
                                    if ((!survivalgene && gol[ij]) || overwrite) {
                                        genecode = golg[nb[k]];
                                        if(coffmod6 > 1) for (genecode1=0ull,k1=0;k1<5;k1++) genecode1 |= ( (genecode>>(40+(k1<<2)+coffmod6-2)) & 0x1ull ) << (k1<<3);
                                        else genecode1 = (genecode >> coffmod6) & 0x0101010101; // prepare lookup in all 5 modules on gene : in this case coff subset bits in right place
                                        genecode &= 0xfcfcfcfcfc; //
                                        genecode |= genecode1;                 // replace lowest 2 bits of 8-bit part of 12-bit modules with appropriate coff subset bits
                                    }
                                    if (!survivalgene && gol[ij]) {
                                        PATTERN8(genecode, pat, found);        //final decision for survival?  NB all 0 seq encodes survival with 0 live nbs
                                        if(repscheme&R_1_nb_OR_AND)
                                            survivenbs |= found? 1ull : 0ull;
                                        else
                                            survivenbs &= found? 1ull : 0ull;
                                    }
                                    if (overwrite) {
                                        PATTERN8(genecode, (0x80|pat), found); //final decision for birth?
                                        if(repscheme&R_1_nb_OR_AND)
                                            birthnbs |= found? 1ull : 0ull;
                                        else
                                            birthnbs &= found? 1ull : 0ull;
                                    }
                                }
                            }
                            if (survivalgene && gol[ij]) {                     // survival determined by central gene in this case
                                genecode = golg[ij];
                                if(coffmod6 > 1) for (genecode1=0ull,k1=0;k1<5;k1++) genecode1 |= ( (genecode>>(40+(k1<<2)+coffmod6-2)) & 0x1ull ) << (k1<<3);
                                else genecode1 = (genecode >> coffmod6) & 0x0101010101; // prepare lookup in all 5 modules on gene : in this case coff subset bits in right place
                                genecode &= 0xfcfcfcfcfc; //
                                genecode |= genecode1;                         // replace lowest 2 bits of 8-bit part of 12-bit modules with appropriate coff subset bits
                                PATTERN8(genecode, pat, found);                //final decision for survival?  NB all 0 seq encodes survival with 0 live nbs
                                survivenbs = found? 1ull : 0ull;
                            }
                            birth &= birthnbs;
                            survive &= survivenbs;
                        }
                        else {                                                 // selection == 14
                            if(repscheme&R_1_nb_OR_AND)
                                for (genecode=0ull,k=0;k<8;k++)                // decodes genes with fixed length encoding by OR
                                    genecode |= (gol[nb[k]]?golg[nb[k]]:0ull); // OR of live neighbours encodes birth rule & survival rule
                            else
                                for (genecode=~0ull,k=0;k<8;k++)               // decodes genes with fixed length encoding by AND
                                    genecode &= (gol[nb[k]]?golg[nb[k]]:~0ull);// AND of live neighbours encodes birth rule & survival rule
                            if(survivalgene) genecode = (genecode&(0xffffffffull<<32)) | (golg[ij]&0xffffffffull); // if central gene determines survival
                            genecode&=(bmask<<32)|smask;                       // actually no longer needed since test done above
                            if (gol[ij]) survive &= (genecode>>(sumoffs[s]+coff))&0x1ull;        // 28.7.2019 changed from |=
                            if (overwrite) birth &= (genecode>>(32+sumoffs[s]+coff))&0x1ull;
                        }
                    }
                }
            }
            else if (rulemodij==2){                                             // hard death on membrane defined above via macro "membrane"
                survive = 0ull;
                birth = 0ull;
            }
            else {
                survive = s2or3;
                birth = s2or3&s&0x1&~gol[ij];
            }
        }
        else survive = birth = s2or3 = gols = 0ull;

        finish_update_ij(ij,s,golij,gols,nb1i,nbmask,nb,survive,birth,gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);
    }  // end for ij

    finish_update(newgol, newgolg, newgolgstats, newgolb, newgolr, nbshist);
}
//------------------------------------------------------------------ genelife_update --------------------------------------------------------------------
void genelife_update (int nsteps, int nhist, int nstat) {
    /* update GoL and gene arrays for toroidal field which has side length which is a binary power of 2 */
    /* encode as much as possible without if structures (use ? : instead) in update routines for optimal vector treatment */
    int k,t,npop;
    uint64_t *newgol, *newgolg, *newgolgstats, *newgolb, *newgolr;
    genedata *genedatap = NULL;
    int totalpoptrace(uint64_t gol[]);                                        // calculate total current population and store in scrolling trace npopulation
    int activitieshash(void);                                                 // count activities of all currently active gene species
    int activitieshashquad(void);                                             // count activities of all currently active quad pattern species
    int get_genealogies(genedata genealogydat[], int narraysize);             // genealogies of all currently active species
    int clonealogies(void);                                                   // clonealogies of all currently active clones
    void tracestats(uint64_t gol[],uint64_t golg[], uint64_t golgstats[], int NN2); // trace statistics based on gol,golg
    void countconfigs(void);
    void countspecies1(uint64_t gol[], uint64_t golg[], int N2);              // count species

    nhistG = nhist;                                                           // intervals for collecting histograms
    nstatG = nstat;
    
    for (t=0; t<nsteps; t++) {                                                // main iteration loop for nsteps
        newgol = planes[newPlane];
        newgolg = planesg[newPlane];
        newgolgstats = planesgs[newPlane];
        newgolb = planesb[newPlane];
        newgolr = planesr[newPlane];
        
        totsteps++;                                                           // simulation step counter
        totdisp++;                                                            // currently every step is counted for display in activities

        if (selection<8)        update_23(gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);           // calculate next iteration with detailed variants of version s=2-3
        else if (selection<10)  update_lut_sum(gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);      // calculate next iteration for lut sum (gene coded)   version s=1-8
        else if (selection<12)  update_lut_dist(gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);     // calculate next iteration for lut dist (corner/edge) version s=1-7
        else if (selection<14)  update_lut_canon_rot(gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);// calculate next iteration for lut canonical rotation version s=2-6
        else if (selection<16)  update_lut_2D_sym(gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);   // calculate next iteration for lut fully 2D symmetric version s=0-4

        if((diagnostics & diag_offset_statistics) && nhist && (totsteps%nhist == 0)) countconfigs(); // count configurations
        if((diagnostics & diag_general_statistics) && nstat && (totsteps%nstat == 0)) tracestats(gol,golg,golgstats,N2); // time trace point

        curPlane = (curPlane+1) % numPlane;                                   // update plane pointers to next cyclic position
        newPlane = (newPlane+1) % numPlane;
        gol = planes[curPlane];                                               // get planes of gol,golg,golb,golr,golgstats data
        golg = planesg[curPlane];
        golgstats = planesgs[curPlane];
        golb = planesb[curPlane];
        golr = planesr[curPlane];
        
        if (diagnostics & diag_scrolling_trace) npop= totalpoptrace(gol);     // calculate total current population and store in scrolling trace npopulation
        
        if ((diagnostics & diag_activities) && (diagnostics & diag_hash_genes)) {
            nspeciesgene=activitieshash();                                    // colors acttrace and sets current population arrays, need to run always for continuity
            if(nspeciesgene<0) fprintf(stderr,"error returned from activitieshash\n");
        }
        if ((diagnostics & diag_activities ) && (diagnostics & diag_hash_patterns)) {
            nspeciesquad=activitieshashquad();                                 // colors acttraceq and sets current population arrays, need to run always for continuity
            if(nspeciesquad<0) fprintf(stderr,"error returned from activitieshashquad\n");
        }
        if(colorfunction==6 || colorfunction==7 || colorfunction2==6 || colorfunction2==7) { // genealogies
            ngenealogydeep=get_genealogies(genedatap,0);                       // calculates and colors genealogytrace
            if(ngenealogydeep<0) fprintf(stderr,"error returned from genealogies\n");
            nclonealogydeep=clonealogies();                                    // calculates and colors clonealogytrace
            if(nclonealogydeep<0) fprintf(stderr,"error returned from clonealogies\n");
        }
        if(!(totsteps%10) && colorupdate1) {
            if (diagnostics & diag_hash_genes)    nallspecies     = hashtable_count(&genetable);
            if (diagnostics & diag_hash_patterns) nallspeciesquad = hashtable_count(&quadtable);
            if (diagnostics & diag_hash_clones)   nallclones = hashtable_count(&clonetable);
            fprintf(stderr,"step %6d:",totsteps);
            if(diagnostics & diag_hash_genes) {
                countspecies1(gol, golg, N2);
                fprintf(stderr," genes %d/%d (extant/all)",nspeciesgene,nallspecies);
            }
            if(diagnostics & diag_hash_patterns) {
                fprintf(stderr," patterns %d/%d (extant/all)",nspeciesquad+nspeciessmall,nallspeciesquad+nallspeciessmall);
            }
            if (diagnostics & diag_hash_clones) {
                fprintf(stderr," clones %d (all)",nallclones);
            }
            fprintf(stderr," ambiguous state s=4 frequency %d",ambigsum);
            fprintf(stderr,"\n");
            fprintf(stderr," neighbour occupation statistics at birth (NW,N,NE,E,SE,S,SW,W):"); for(k=0;k<8;k++) fprintf(stderr," %d",nbshist[k]); fprintf(stderr,"\n");
            fprintf(stderr,"__________________________________________________________________________________________________________________________________________\n");
        }

    }
}
//----------------------------------------------------------------- initialize_planes -------------------------------------------------------------------
void initialize_planes(int offs[],  int Noffsets) {
    int i,j,idx;
    static int notfirst = 0;

    curPlane = 0;
    newPlane = 1;
    if (notfirst)   return;     // need to fix memory free at two levels unless this fix: no changes in planes structure during run allowed
    notfirst = 1;

    planes[0]  = plane0;  planes[1]  = plane1;     // initialize plane pointers:
    planesg[0] = planeg0; planesg[1] = planeg1;
    planesgs[0] = planegs0; planesgs[1] = planegs1;
    planesb[0]  = planeb0;  planesb[1]  = planeb1;
    planesr[0]  = planer0;  planesr[1]  = planer1;
#if maxPlane > 2
    planes[2]  = plane2;  planes[3]  = plane3;
    planesg[2] = planeg2; planesg[3] = planeg3;
    planesgs[2] = planegs2; planesgs[3] = planegs3;
    planesb[2]  = planeb2;  planesb[3]  = planeb3;
    planesr[2]  = planer2;  planesr[3]  = planer3;
#endif
#if maxPlane > 4
    planes[4]  = plane4;  planes[5]  = plane5;  planes[6]  = plane6;  planes[7]  = plane7;
    planesg[4] = planeg4; planesg[5] = planeg5; planesg[6] = planeg6; planesg[7] = planeg7;
    planesgs[4] = planegs4; planesgs[5] = planegs5; planesgs[6] = planegs6; planesgs[7] = planegs7;
    planesb[4]  = planeb4;  planesb[5]  = planeb5;  planesb[6]  = planeb6;  planesb[7]  = planeb7;
    planesr[4]  = planer4;  planesr[5]  = planer5;  planesr[6]  = planer6;  planesr[7]  = planer7;
#endif

    if (!(diagnostics & diag_offset_statistics)) return;
    
    if(Noffsets%3 !=0) fprintf(stderr,"Size of offsets array not a multiple of 3 as expected.");
    Noff = Noffsets/3;        // Noff global
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
    int toff, tmn = Noffsets;
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
}
//---------------------------------------------------------------- readFile and writeFile ---------------------------------------------------------------
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
int writeFile(char *fileName)  {            // initialize 32x32 genepat file with all empty sites
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
//---------------------------------------------------------------- initialize ---------------------------------------------------------------------------
void testmacros() {
    int j,k;
    uint64_t g;
    g=1ull;                                                              // test of FIRST1INDEX
    for(j=0;j<10;j++) {
        FIRST1INDEX(g,k);
        fprintf(stderr,"test of first1index cnt %d val %llx index %d\n",j,g,k);
        g*=42;
    }

    g=1ull;                                                              // test of PATTERN4
    for(k=0;k<10;k++) {
        int found;
        for(j=0;j<16;j++) {
            PATTERN4(g,j,found);
            fprintf(stderr,"test of pattern4 pat %x val %llx found? %d\n",j,g,found);
        }
        g*=42;
    }

    g=1ull;                                                              // test of PATTERN8
    for(k=0;k<10;k++) {
        int found;
        for(j=0;j<64;j++) {
            PATTERN8(g,j,found);
            fprintf(stderr,"test of pattern8 pat %x val %llx found? %d\n",j,g,found);
        }
        g*=42;
    }
    
    for (j=0;j<257;j++) {
        fprintf(stderr,"test of patterns j %5d log2upper(j) %5d log2upper(sqrtupper(j)) %5d\n",j,log2upper(j),log2upper(sqrtupper(j)));
    }
}
//........................................................................................................................................................
void initialize(int runparams[], int nrunparams, int simparams[], int nsimparams) {
    int hcnt;
    int ij,ij1,i0,j0,i,j,Nf,k,cnt,icf,nstartgenes;
    unsigned int ncodingin;
    uint64_t g;
    static unsigned int rmask = (1 << 15) - 1;
    static int notfirst = 0;
    uint64_t startgenes[16];
    char *golgin;

    rulemod = runparams[0];
    repscheme = runparams[1];
    selection = runparams[2];
    overwritemask = runparams[3];
    survivalmask = runparams[4];
    colorfunction = runparams[5];
    initfield = runparams[6];
    birthmask=runparams[7];
    ancselectmask=runparams[8];
    colorfunction2 = runparams[9];

    pmutmask = (unsigned int) simparams[0];                                      // low values of pmutmask <32 are interpreted as -nlog2pmut
    if(pmutmask<32) pmutmask = (0x1 << (32-pmutmask)) - 0x1ull;                  // NB if pmut==0, pmutmask==zero, no mutation.
    initial1density = simparams[1];
    initialrdensity = simparams[2];
    
    ncodingin = simparams[3];                                                    // used in selection  4,5,6,8,
    ncoding = ncodingin & 0xff;
    ncoding2 = (ncodingin>>8) & 0xff;
    if (ncoding > 64) { fprintf(stderr,"value %d of ncoding is out of range\n",ncoding);ncoding = 64;}
    if (selection==8 || selection == 9) if (ncoding<1 || ncoding>4) ncoding = 4; // ncoding range restriction for selection = 8 ie 16*ncoding bits of gene used
    if (selection<8) codingmask = (1ull<<ncoding2)-1ull;                         // bit mask corresponding to ncoding2 bits, only used in connection with add2ndmask1st
    else if (selection<10) codingmask = (1ull<<ncoding)-1ull;                    // coding mask used to encode number of bits per LUT (1-4)
    
    startgenechoice = simparams[4];
    if(nsimparams > 5) ranseed = simparams[5];
    
    srand(ranseed); // Range: rand returns numbers in the range of [0, RAND_MAX ), and RAND_MAX is specified with a minimum value of 32,767. i.e. 15 bit
    randstate[0] = rand();randstate[1] = rand();                                 // state vector for dedicated 64-bit random number generator macro RAND128P
    cnt = 0;
    totsteps = 0;
    totdisp = 0;
    statcnts = 0;
    quadrants = -1;
    rbackground = 0;
    quadcollisions = 0;
    randominflux = 0;
    vscrolling = last_scrolled = vscrolly = 0;
    quadrants = -1;
    gene0=0ull;                                                                 // normally default gene is 0ull : unused when gol state not live
    nstartgenes = 8;
    
    // testmacros();                                                             // test macros used to accelerate processing
    // writeFile("genepat.dat");                                                 // can be used to initialize formatted template for gene input of 32x32 array

    fprintf(stderr,"___________________________________________________________________________________________\n");
    fprintf(stderr,"_________________________________ genelife simulation _____________________________________\n");
    fprintf(stderr,"runparams %d %d %d %d %d %d %d %d %d %d\n",runparams[0],runparams[1],runparams[2],
                    runparams[3],runparams[4],runparams[5],runparams[6],runparams[7],runparams[8],runparams[9]);
    fprintf(stderr,"simparams %d %d %d %d %d %d\n",simparams[0],simparams[1],simparams[2],simparams[3],simparams[4],ranseed);
    fprintf(stderr,"pmutmask %x (NB 0 means no mutation)\n",pmutmask);

    switch (selection) {                                                         // initialize starting genes depending on selection model, encoding and symmetry
        case 0:  for (k=0;k<4;k++) { startgenes[k] = 0xf0f0f0f0f0f0f0f0; startgenes[k+4] = 0x0f0f0f0f0f0f0f0f;} break;
        case 1:  for (k=0;k<8;k++)   startgenes[k] = ((0x1ull<<k*3)-1ull)<<20;break;
        case 2:  for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | 0x606ull;break; // to DEBUG 2-8 difference otehrwise same as case 3:
        case 3:  for (k=0;k<8;k++)   startgenes[k] =(((0x1ull<<20)-1ull)<<20)+((0x1ull<<k)-0x1ull);break;
        case 4:  for (k=0;k<8;k++)  {g = 0xff0ull; startgenes[k] = k<4 ? g+1 : ~((g<<16));} break;
        case 5:  for (k=0;k<8;k++)  {g = 0xf0ull + k; startgenes[k] = k<4 ? g : (~g)|(0xfull<<16);} break;
        case 6:
        case 7:  for (k=0;k<8;k++) startgenes[k]=(0x1ull<<(4+k*8))-1ull; break;

        case 8:  if (((repscheme>>4)&0x7)==7)
                      genegol[selection-8] = (codingmask<<((8+3-1)*ncoding))|(codingmask<<((8+2-1)*ncoding))|(codingmask<<((2-1)*ncoding))|(codingmask<<((3-1)*ncoding));
                 else genegol[selection-8] = (codingmask<<((8+3-1)*ncoding))|(codingmask<<((2-1)*ncoding))|(codingmask<<((3-1)*ncoding));
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];break;    // put up to 3 extra bits at top to ensure all nr 1s values occupied
        case 9:  if (((repscheme>>4)&0x7)==7) genegol[selection-8] = 0xba32ull;                             // extended rule for scissors-paper-stone-well gliders
                 else                         genegol[selection-8] = 0xb32ull;                              //GoL rule for survival in totalistic LUT case, variable length encoding
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];break;    // put up to 3 extra bits at top to ensure all nr 1s values occupied
        case 10: if (((repscheme>>4)&0x7)==7) genegol[selection-8] = (0x7ull<<2)|(0xfull<<5)|(0x7ull<<34)|(0xfull<<37); //GoL rule for S23 B23 in corner/edge dist LUT case, fixed length encoding
                 else                         genegol[selection-8] = (0x7ull<<2)|(0xfull<<5)|(0xfull<<37);  //GoL rule for S23 B3 in corner/edge dist LUT case, fixed length encoding
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];break;
        case 11: if (((repscheme>>4)&0x7)==7) genegol[selection-8] = 0xbfa73f27ull;                         //GoL rule for S23 B23 in corner/edge dist LUT case, variable length encoding
                 else                         genegol[selection-8] = 0xbf3f27ull;                           //GoL rule for S23 B3 in corner/edge dist LUT case, variable length encoding
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];break;
        case 12: if (((repscheme>>4)&0x7)==7) genegol[selection-8] = 0xfull|(0x7full<<4)|(0xfull<<32)|(0x7full<<36); //GoL rule for S23 B23 in canonical rotation case, fixed length encoding
                 else                         genegol[selection-8] = 0xfull|(0x7full<<4)|(0x7full<<36);     //GoL rule for 2,3 s and 3 b in canonical rotation case:  4,7,10,7,4 configs s=2,3,4,5,6
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];break;
        case 13: if (((repscheme>>4)&0x7)==7) genegol[selection-8] = 0x04cda7bbb73b3727;                    //GoL rule for S23 B23 with modular encoding for canonical rotation case
                 else                         genegol[selection-8] = 0x00cd00bbb73b3727;                    //GoL rule for S23 B3 with modular encoding for canon rot'n case: 4,7,10,7,4 configs s=2,3,4,5,6
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];break;
        case 14: if (((repscheme>>4)&0x7)==7) genegol[selection-8] = (0x3full<<3)|(0x3ffull<<9)|(0x3full<<(32+3))|(0x3ffull<<(32+9)); //GoL rule for S23 B23 for 2D_sym case, fixed length encoding
                 else                         genegol[selection-8] = (0x3full<<3)|(0x3ffull<<9)|(0x3ffull<<(32+9)); //GoL rule for S23 B3 for 2D_sym case with 1,2,6,10,13 configs for s=0,1,2,3,4
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];break;
        case 15: if (((repscheme>>4)&0x7)==7) genegol[selection-8] = 0x03ffffb7b3a33323;                    //GoL rule for S2(all)3(1st 6) B23 with modular encoding for 2D_sym case (NB smask 0 anyway)
                 else                         genegol[selection-8] = 0x03f3ffb7b3373323;                    //GoL rule for S23 B3 with modular encoding for 2D_sym case for s=0,1,2,3,4
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];break;

        default: for (k=0;k<8;k++) startgenes[k]=(0x1ull<<(4+k*8))-1ull;
    }
    
    if (diagnostics & diag_general_statistics) {                                  // general statistics
        if (livesites !=NULL) {free(livesites);livesites = NULL;}
        if (genestats !=NULL) {free(genestats);genestats = NULL;}
        if (stepstats !=NULL) {free(stepstats);stepstats = NULL;}
        if (configstats != NULL) {free(configstats);configstats = NULL;}
    
        arraysize = startarraysize;
        livesites = (int *) malloc(arraysize * sizeof(int));
        genestats = (int *) malloc(arraysize * 4 * sizeof(int));
        stepstats = (int *) malloc(arraysize * 10 * sizeof(int));
        if (nhistG==nstatG) configstats = (int *) malloc(arraysize * Noff * sizeof(int));
    }
    
    curPlane = 0;                                                                 // if we rerun initialize, we want to restart plane cycling from zero
    newPlane = 1;
    gol = planes[curPlane];
    golg = planesg[curPlane];
    golgstats = planesgs[curPlane];
    golb = planesb[curPlane];
    golr = planesr[curPlane];
    
    if(notfirst) {
        if(diagnostics & diag_hash_genes) hashtable_term(&genetable);
        if(diagnostics & diag_hash_patterns) {
            hashtable_term(&quadtable);
            memset(smallpatts,0,sizeof(smallpatt)*65536);
            // for (ij=0;ij<65536;ij++) smallpatts[ij].activity = 0;              // initialize small pattern table to no patterns hit, already done
        }
        if(diagnostics & diag_hash_clones) hashtable_term(&clonetable);
    }
    
    if(diagnostics & diag_hash_genes) hashtable_init(&genetable,sizeof(genedata),N2<<2,0);   // initialize dictionary for genes
    if(diagnostics & diag_hash_patterns) hashtable_init(&quadtable,sizeof(quadnode),N2<<2,0);// initialize dictionary for quadtree patterns
    if(diagnostics & diag_hash_clones) hashtable_init(&clonetable,sizeof(clonedata),N2<<4,0);// initialize dictionary for clones
    
    notfirst = 1;
    if (initfield==1) {                              // input from file genepat.dat with max size of 32*32 characters
        golgin = (char *) malloc(32* 32 * sizeof(char));
        icf=readFile(golgin, "genepat.dat");
        if (icf != 32*32) {
            icf = 0;
            fprintf(stderr,"error reading file, %d not 32*32 chars\n",icf);
        }
        for (ij=0; ij<N2; ij++) {
            gol[ij] = 0ull;
            golg[ij] = gene0;
            golgstats[ij] = 0ull;
            golb[ij] = 0ull;
            golr[ij] = 0ull;
        }
        for (ij1=0; ij1<32*32; ij1++) {
            if(N<32) {fprintf(stderr,"Error, array dim %d too small for file array dim %d\n",N,32);break;}
            ij=(N>>1)-16+(ij1&0x1f)+ N*((N>>1)-16+(ij1>>5));
            if (golgin[ij1] > 0)    {                // if live cell
                gol[ij] = 1ull;
                if(golgin[ij1] <= 8 ) golg[ij] = startgenes[golgin[ij1]-1];
                else if (golgin[ij1]>='0' && golgin[ij1]<'8') golg[ij] = startgenes[golgin[ij1]-'0'];
                else golg[ij] = startgenes[7];
                cnt++;
                golb[ij] = rootclone + ij;               // initialize clone to new clone at time 0
            }
            // if (golg[ij] == 0 && gol[ij] != 0) fprintf(stderr,"zero gene at %d\n",ij);
        }

    }
    else if (initfield>=0) {                         // initfield gives linear size of random block for initialization (0 => full frame, as before)
        Nf = initfield;
        if (Nf==0 || Nf>N) Nf=N;
        for (ij=0; ij<N2; ij++) {
            gol[ij] = 0ull;
            golg[ij] = gene0;
            golgstats[ij] = 0ull;
            golb[ij] = 0ull;
            golr[ij] = 0ull;
        }
        i0 = j0 = (N>>1)-(Nf>>1);
        for (i=0; i<Nf; i++) {
            for (j=0; j<Nf; j++) {
                ij=i0+i+N*(j0+j);
                gol[ij] |= ((rand() & rmask) < initial1density)?1ull:0ull;
            }
        }
        for (ij=0; ij<N2; ij++) {
            g = 0ull;
            if (gol[ij]) {                          //  fill with random genome g or randomly chosen startgene depending on initialrdensity
                if (((unsigned) rand() & rmask) < initialrdensity) {for (k=0; k<64; k++) g = (g << 1) | (rand() & 0x1);g=gene0^g;}
                else if (startgenechoice == nstartgenes) g = startgenes[0xf & rand() & (nstartgenes-1)];
                else if (startgenechoice > nstartgenes) fprintf(stderr,"startgenechoice %d out of range\n",startgenechoice);
                else g = startgenes[0xf & startgenechoice & (nstartgenes-1)];
                cnt++;
            }
            golg[ij] = g;
            golb[ij] = rootclone + ij;
        }
    }
    else {                                          // initfield < 0, use array values from stashed
        for (ij=0; ij<N2; ij++) {
            gol[ij] = stashgol[ij];
            golg[ij] = stashgolg[ij];
            golgstats[ij] = 0ull;
            golb[ij] = rootclone + ij;              // might want to define stashgolb and use it here
            golr[ij] = 0ull;
        }
    }

    if(diagnostics & diag_scrolling_trace) {
        for (i=0;i<N;i++) npopulation[i]=0;         // initialize scrolling total population trace
        for (ij=0; ij<N2; ij++) {
            poptrace[ij]=rootgene;                  // initialize population traces to root gene
            acttrace[ij]=rootgene;                  // initialize activity traces to root gene
            acttraceq[ij]=rootgene;                 // initialize activity traces of patterns to root gene
            acttraceqt[ij]=0;                       // initialize activity traces of types of patterns to 0
            genealogytrace[ij] = rootgene;          // initialize genealogy traces to root gene
        }
    }

    if(diagnostics & diag_longtime_trace) {
        for (i=0;i<N*nNhist;i++) npopulation1[i]=0; // initialize longer beginning total population trace
        for (ij=0; ij<N2*nNhist; ij++) {
            poptrace1[ij]=rootgene;                 // initialize long population traces to root gene
            acttrace1[ij]=rootgene;                 // initialize long activity traces to root gene
            acttraceq1[ij]=rootgene;                // initialize long activity traces of patterns to root gene
            acttraceqt1[ij]=0;                      // initialize long activity traces of types of patterns to 0
        }
    }
    
    if(diagnostics & diag_hash_genes) {
        for (ij=0; ij<N2; ij++) {
            if(gol[ij]) {
                hashaddgene(ij,golg[ij],rootgene,golb+ij,rootclone+ij,0x1ull);   // totsteps=0 so parentid is spatial ij + flag rootclone
            }
        }
                                                // enumerate gene hash keys and values
        hcnt=hashtable_count(&genetable);
        genotypes = hashtable_keys( &genetable );
        fprintf(stderr,"population size %d with %d different genes\n",cnt,hcnt);
    }

    // if(diagnostics & diag_hash_patterns) qimage = quadimage(newgol,&patt,log2N); // quadtree hash of entire image
    if(diagnostics & diag_component_labels) ncomponents=extract_components(gol);
}
//-------------------------------------------------------------------- set ...---------------------------------------------------------------------------
void set_colorfunction(int colorfunctionin) {
    if((colorfunctionin>12) || (colorfunctionin<0)) fprintf(stderr,"error colorfunction value passed %d out of range\n",colorfunctionin);
    else     colorfunction = colorfunctionin;
}
//.......................................................................................................................................................
void set_colorfunction2(int colorfunctionin) {
    if((colorfunctionin>12) || (colorfunctionin<-1)) fprintf(stderr,"error colorfunction value passed %d out of range\n",colorfunctionin);
    else     colorfunction2 = colorfunctionin;
}
//.......................................................................................................................................................
int setget_act_ymax(int actymax) {                  // sets ymax for activities only if argument nonzero, reads old value
    int ymaxold;
    ymaxold = ymax;
    ymax = actymax;
    return(ymaxold);
}
//.......................................................................................................................................................
int setget_act_ymaxq(int actymaxq) {                  // sets ymax for activities only if argument nonzero, reads old value
    int ymaxqold;
    ymaxqold = ymaxq;
    ymaxq = actymaxq;
    return(ymaxqold);
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
        repscheme |= R_16_quadrant_sele<<quadrant;                          // assumes quadrant selectors are 7 successive bits following R_16_...
    }
}
//.......................................................................................................................................................
void set_randominflux(int randominfluxin) {
    randominflux=randominfluxin;
}
//.......................................................................................................................................................
void set_rbackground(int rbackgroundin, int randominfluxin) {
    rbackground=rbackgroundin;
    randominflux=randominfluxin;
    if(!rbackground) randominflux=0;   // do not leave patch randominflux variable active when turning off random background
}
//.......................................................................................................................................................
unsigned int set_repscheme_bits(int quadrant, int x, int y, unsigned int surviveover[]) {
    unsigned int quadrantval;

    quadrantval=(x<(Nmask>>1)? 0 : 1) + (y<(Nmask>>1)? 0 : 2);              // determine selected quadrant
    if(quadrants >= 0 && quadrants < 5) {                                   // assumes repscheme bits in pairs starting from bit 0 matching quadrants
        repscheme &=  ~(0x3llu<<(quadrants<<1));
        repscheme |=  quadrantval<<(quadrant<<1);
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
//.......................................................................................................................................................
void set_noveltyfilter() {
    noveltyfilter=1-noveltyfilter;
}
//.......................................................................................................................................................
void set_activity_size_colormode() {
    activity_size_colormode = (activity_size_colormode+1) % 4;
}
//.......................................................................................................................................................
void set_gcolors() {
    gcolors = (gcolors+1)%10;
}
//.......................................................................................................................................................
void set_seed(int seed) {
    ranseed = seed;
}
//.......................................................................................................................................................
void set_nbhist(int nbhistin) {
    if(nbhist<nNhist*2) nbhist=nbhistin;
    else fprintf(stderr,"nbhist out of range %d > %d\n",nbhistin,nNhist*2-1);
}
//.......................................................................................................................................................
void set_genealogycoldepth(int genealogycoldepthin) {
    genealogycoldepth = genealogycoldepthin;
}
//.......................................................................................................................................................
void set_ancestortype(int ancestortypein) {
    if(ancestortypein <3) ancestortype = ancestortypein;
    else fprintf(stderr,"ancestor type %d out of range [0..2]\n",ancestortypein);
}
//.......................................................................................................................................................
void set_stash(){               // stash current gol,golg
    int ij;
    for (ij=0; ij<N2; ij++) {
        stashgol[ij] = planes[curPlane][ij];
        stashgolg[ij] = planesg[curPlane][ij];
    }
}
//.......................................................................................................................................................
void set_info_transfer_h(int do_info_transfer, int nbhood) {
    info_transfer_h = do_info_transfer;
    if(nbhood == 3 || nbhood == 5 || nbhood == 7)
        it_nbhood = nbhood;
    else fprintf(stderr,"error in nbhood value %d, allowed values are 3,5,7\n",nbhood);
}
//.......................................................................................................................................................
void set_activityfnlut(int activityfnlutin) {
    activityfnlut = activityfnlutin;
}
//.......................................................................................................................................................
void set_colorupdate1(int update1) {
    colorupdate1 = update1;
}
//------------------------------------------------------------------- get ... ---------------------------------------------------------------------------
int get_log2N() {
    return(log2N);
}
//.......................................................................................................................................................
void get_stash(){               // retrieve current gol,golg from stashed values
    int ij;
    for (ij=0; ij<N2; ij++) {
        planes[curPlane][ij] = stashgol[ij];
        planesg[curPlane][ij] = stashgolg[ij];
    }
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
void get_curgolbr(uint64_t outgolb[], uint64_t outgolr[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolb[ij] = planesb[curPlane][ij];
    }
    for (ij=0; ij<NN; ij++) {
        outgolr[ij] = planesr[curPlane][ij];
    }
}
//.......................................................................................................................................................
void get_stats(int outstats[], int outgtypes[], int outstepstats[], int outconfigstats[], int numStats ){
    int i;
    
    if (!(diagnostics & diag_general_statistics)) {
        fprintf(stderr,"statistics collection not enabled in C\n");
        return;
    }
    if(numStats > arraysize){
        fprintf(stderr,"Ack! numStats = %d  > arraysize = %d\n",numStats,arraysize);
        exit(1);
    }
    for(i=0; i<numStats; i++) outstats[i] = livesites[i];
    for(i=0; i<4*numStats; i++) outgtypes[i] = genestats[i];
    for(i=0; i<10*numStats; i++) outstepstats[i] = stepstats[i];
    if (nhistG==nstatG) for(i=0; i<Noff*numStats; i++) outconfigstats[i] = configstats[i];
}
//.......................................................................................................................................................
void get_acttrace(uint64_t outgolg[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolg[ij] = acttrace[ij];
    }
}
//.......................................................................................................................................................
void get_acttraceq(uint64_t outgolg[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolg[ij] = acttraceq[ij];
    }
}
//.......................................................................................................................................................
void get_poptrace(uint64_t outgolg[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolg[ij] = poptrace[ij];
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
int get_nlive() {
    int ij,nlive;
    for(nlive=0,ij=0;ij<N2;ij++) nlive+= (gol[ij]>0) ? 1 : 0;
    return(nlive);
}
//.......................................................................................................................................................
int get_genealogydepth() {
    int j, jmax, i, nspecies, nspeciesnow;
    uint64_t gene, ancgene;
    int *gindices;
    uint64_t *genes;

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );

    for (i=nspeciesnow=0; i<nspecies; i++)
        if(geneitems[i].popcount) nspeciesnow++;


    gindices = (int *) malloc(nspecies*sizeof(int));
    for (i=j=0; i<nspecies; i++) {
        if(geneitems[i].popcount) {
            gindices[j]=i;
            j++;
        }
        else gindices[nspeciesnow+i-j]=i;
    }

    genes = (uint64_t *) malloc(nspeciesnow*sizeof(uint64_t));

    for (i=0; i<nspeciesnow; i++) {
        genes[i]=genotypes[gindices[i]];
    }
    
    if(ancestortype>0) fprintf(stderr,"Warning: get_genealogydepth currently only implemented for ancestortypes 0 called with %d\n",ancestortype);
    for (i=jmax=0; i<nspeciesnow; i++) {                            // calculate max depth in genealogy jmax
        gene=genes[i];
        ancgene=geneitems[gindices[i]].firstancestor;
        for (j=1;;j++) {
            gene=ancgene;
            if(gene==rootgene) break;                               // reached root, exit j loop
            else {
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
                    ancgene=genedataptr->firstancestor;
                }
                else fprintf(stderr,"ancestor not found in genealogies\n");
            }
        }
        if (j>jmax) jmax=j;
    }
    genealogydepth = jmax+1;
    free(genes); free(gindices);
    return(genealogydepth);
}
//.......................................................................................................................................................
int get_curtime(){
    return(totsteps);
}
//.......................................................................................................................................................
void get_histo(int outhisto[],int numHistoC){
    int i;
   
    if (!(diagnostics & diag_offset_statistics)) {
        fprintf(stderr,"histogram of offsets not enabled, activate diag_offset_statistics\n");
        return;
    }
    if(numHistoC != numHisto){
        fprintf(stderr,"Ack! numHisto = %d  != numHistoC = %d\n",numHisto,numHistoC);
        exit(1);
    }
    for(i=0; i<numHisto; i++) outhisto[i] = histo[i];
}
//.......................................................................................................................................................
int get_activities(uint64_t actgenes[], int activities[], int narraysize) {
    int k, nlivegenes, nspecies;

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );
    // fprintf(stderr,"The number of different species that have ever existed is %d\n",nspecies);

    for (k=nlivegenes=0; k<nspecies; k++) {
        if((genedataptr = (genedata *) hashtable_find(&genetable, genotypes[k])) != NULL) {
            if(genedataptr->popcount) {
                if (nlivegenes <= narraysize) {
                    actgenes[nlivegenes] = genotypes[k];
                    activities[nlivegenes] = genedataptr->activity;
                }
                nlivegenes++;
            }
        }
        else fprintf(stderr,"get_activities error, no entry for gene %llx in hash table\n", genotypes[k]);
    }
    if (nlivegenes > narraysize) fprintf(stderr,"Error: array size %d to small to hold live activities %d, increase it\n",narraysize,nlivegenes);

    return nlivegenes;
}
//.......................................................................................................................................................
int get_all_activities(uint64_t genes[], int activities[], int narraysize) {
// get_all_activities   get all activity statistics of genes (since t=0) from C to python
    int k, nspecies;

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata *) hashtable_items( &genetable );
    // fprintf(stderr,"The number of different species that have ever existed is %d\n",nspecies);
    if (nspecies > narraysize) {
        fprintf(stderr,"Error: array size %d to small to hold all activities %d, increase it\n",narraysize,nspecies);
        return nspecies;
    }

    for (k=0; k<nspecies; k++) {
        if((genedataptr = (genedata *) hashtable_find(&genetable, genotypes[k])) != NULL) {
            genes[k] = genotypes[k];
            activities[k] = genedataptr->activity;
        }
        else fprintf(stderr,"get_all_activities error, no entry for gene %llx in hash table\n", genotypes[k]);
    }
    return nspecies;
}
//.......................................................................................................................................................
int get_quad_activities(uint64_t quads[], int activities[], int narraysize) {
// get_quad_activities  get *live* activity statistics of quads (since t=0) from C to python
    int k, nspecies, livecnt;
    quadnode *q;

    nspecies = hashtable_count(&quadtable);
    quadkeys = hashtable_keys(&quadtable);
    quaditems = (quadnode *) hashtable_items( &genetable );

    for (k=0,livecnt = 0; k<nspecies; k++)
        if((q = (quadnode *) hashtable_find(&quadtable, quadkeys[k])) != NULL){
            if(q->lasttime == totsteps) // test for currently live
                livecnt++;
        } else {
            fprintf(stderr,"get_quad_activities error, no entry for quad %llx in hash table\n", quadkeys[k]);
        }
    // fprintf(stderr,"The number of different species that have ever existed is %d\n",nspecies);
    if (livecnt > narraysize) {
        fprintf(stderr,"Error: array size %d to small to hold all quad activities %d, increase it\n",narraysize,livecnt);
        return livecnt;
    }

    for (k=0; k<nspecies; k++) {
        if((q = (quadnode *) hashtable_find(&quadtable, quadkeys[k])) != NULL) {
            if(q->lasttime == totsteps){ // test for currently live
                quads[k] = quadkeys[k];
                activities[k] = q->activity;
            }
        }
        else fprintf(stderr,"get_quad_activities error, no entry for quad %llx in hash table\n", quadkeys[k]);
    }
    return nspecies;
}
//.......................................................................................................................................................
int get_all_quad_activities(uint64_t quads[], int activities[], int narraysize) {
// get_all_quad_activities  get all activity statistics of quads (since t=0) from C to python
    int k, nspecies;
    quadnode *q;

    nspecies = hashtable_count(&quadtable);
    quadkeys = hashtable_keys(&quadtable);
    quaditems = (quadnode *) hashtable_items( &genetable );
    // fprintf(stderr,"The number of different species that have ever existed is %d\n",nspecies);
    if (nspecies > narraysize) {
        fprintf(stderr,"Error: array size %d to small to hold all quad activities %d, increase it\n",narraysize,nspecies);
        return nspecies;
    }

    for (k=0; k<nspecies; k++) {
        if((q = (quadnode *) hashtable_find(&quadtable, quadkeys[k])) != NULL) {
            quads[k] = quadkeys[k];
            activities[k] = q->activity;
        }
        else fprintf(stderr,"get_quad_activities error, no entry for quad %llx in hash table\n", quadkeys[k]);
    }
    return nspecies;
}
//.......................................................................................................................................................
int get_small_activities(uint64_t smalls[], int activities[], int narraysize) {
// get_small_activities  get *live only* activity statistics of smallpatts (since t=0) from C to python
    int k, nspecies;

    if (narraysize<65536) {
        fprintf(stderr,"Error in get_small_activities : called with insufficent smallpatt holding array size %d < %d\n",narraysize,65536);
        return -1;
    }
    nspecies =0;
    for (k=0; k<65536; k++) {
        if(smallpatts[k].lasttime == totsteps){
            nspecies += smallpatts[k].activity ? 1 : 0;
            activities[k] = smallpatts[k].activity;
        }
    }
    return nspecies;
}
//.......................................................................................................................................................
int get_all_small_activities(uint64_t smalls[], int activities[], int narraysize) {
// get_all_small_activities  get all activity statistics of quads (since t=0) from C to python
    int k, nspecies;

    if (narraysize<65536) {
        fprintf(stderr,"Error in get_small_activities : called with insufficent smallpatt holding array size %d < %d\n",narraysize,65536);
        return -1;
    }
    nspecies =0;
    for (k=0; k<65536; k++) {
        nspecies += smallpatts[k].activity ? 1 : 0;
        activities[k] = smallpatts[k].activity;
    }
    return nspecies;
}
//.......................................................................................................................................................
int get_connected_comps(unsigned int outlabel[], unsigned int outconnlen[], int x, int y) {
    int i,ij;
    for (ij=0; ij<N2; ij++) {
        outlabel[ij] = (unsigned int) label[ij];
    }
    for (i=1;i<ncomponents+1;i++) {
        outconnlen[i] = (unsigned int) connlen[i];
    }
    xdisplay = x;
    ydisplay = y;
    return ncomponents;
}
//.......................................................................................................................................................
int get_ncomponents() {
    return(ncomponents);
}
//.......................................................................................................................................................
int get_components(component components[],int narraysize) {
    int i;
    if (narraysize<ncomponents) {
        fprintf(stderr,"Error in get_components : called with insufficent component holding array size %d < %d\n",narraysize,ncomponents);
        return -1;
    }
    for (i=1;i<=ncomponents;i++) {
        components[i-1]=complist[i];
        /*components[i-1].N=complist[i].N;
        components[i-1].S=complist[i].S;
        components[i-1].W=complist[i].W;
        components[i-1].E=complist[i].E;
        components[i-1].lastrc=complist[i].lastrc;
        components[i-1].label=complist[i].label;
        components[i-1].log2n=complist[i].log2n;
        components[i-1].patt=complist[i].patt;
        components[i-1].quad=complist[i].quad;
        components[i-1].pixels=complist[i].pixels;
        components[i-1].gcolor=complist[i].gcolor;*/
        /* if (i<100) fprintf(stderr,"Component %d: (N,W,S,E)=(%d,%d,%d,%d), lastrc=%d, label=%d, log2n=%d, patt=%d, quad=%llx, pixels=%d\n",i,
            complist[i].N,complist[i].W,complist[i].S,complist[i].E,complist[i].lastrc,complist[i].label,
            complist[i].log2n,complist[i].patt,complist[i].quad,complist[i].pixels); */
    }

    return ncomponents;
}
//.......................................................................................................................................................
int get_smallpatts(smallpatt smallpattsout[],int narraysize) {
    int i,count;
    if (narraysize<65536) {
        fprintf(stderr,"Error in get_smallpatts : called with insufficent smallpatt holding array size %d < %d\n",narraysize,65536);
        return -1;
    }
    for (count=i=0;i<65536;i++) {
        // smallpattsout[i].size= smallpatts[i].size;
        smallpattsout[i].topactivity= smallpatts[i].topactivity;
        smallpattsout[i].activity= smallpatts[i].activity;
        count += smallpatts[i].activity ? 1 : 0;
        smallpattsout[i].firsttime= smallpatts[i].firsttime;
        smallpattsout[i].lasttime= smallpatts[i].lasttime;
    }
    return count;
}
//.......................................................................................................................................................
int get_quadnodes(quadnode quadnodes[],int narraysize) {
    int i;

    // these three calls executed through hashactivityquad if colorfunction 9 or 10
    // nallspeciesquad = hashtable_count(&quadtable);
    // quadkeys = hashtable_keys(&quadtable);
    // quaditems = (quadnode*) hashtable_items( &quadtable );

    if (narraysize<nallspeciesquad) {
        fprintf(stderr,"Error in get_quadnodes : called with insufficent quadnode holding array size %d < %d\n",narraysize,nallspeciesquad);
        return -1;
    }

    for (i=0;i<nallspeciesquad;i++) {
        quadnodes[i]=quaditems[i];
        /*quadnodes[i].hashkey=quaditems[i].hashkey;
        quadnodes[i].nw=quaditems[i].nw;
        quadnodes[i].ne=quaditems[i].ne;
        quadnodes[i].sw=quaditems[i].sw;
        quadnodes[i].se=quaditems[i].se;
        quadnodes[i].isnode=quaditems[i].isnode;
        quadnodes[i].size=quaditems[i].size;
        quadnodes[i].activity=quaditems[i].activity;
        quadnodes[i].pop1s=quaditems[i].pop1s;
        quadnodes[i].firsttime=quaditems[i].firsttime;
        quadnodes[i].lasttime=quaditems[i].lasttime;
        quadnodes[i].topactivity=quaditems[i].topactivity;*/
    }

    return nallspeciesquad;
}
//.......................................................................................................................................................
int get_genes(genedata genelist[],int narraysize) {
    int i;

    // these three calls executed already through hashactivity
    // nallspecies = hashtable_count(&genetable);
    // genotypes = hashtable_keys(&genetable);
    // geneitems = (genedata*) hashtable_items( &genetable );

    if (narraysize<nallspecies) {
        fprintf(stderr,"Error in get_genes : called with insufficent genedata holding array size %d < %d\n",narraysize,nallspecies);
        return -1;
    }

    for (i=0;i<nallspecies;i++) {
        genelist[i]=geneitems[i];          /* shallow copy OK if no pointers being copied, otherwise the structures they point to will not be copied */
        /* genelist[i].popcount=geneitems[i].popcount;
        genelist[i].firsttime=geneitems[i].firsttime;
        genelist[i].lasttime=geneitems[i].lasttime;
        genelist[i].lastextinctiontime=geneitems[i].lastextinctiontime;
        genelist[i].activity=geneitems[i].activity;
        genelist[i].nextinctions=geneitems[i].nextinctions;
        genelist[i].gene=geneitems[i].gene;
        genelist[i].firstancestor=geneitems[i].firstancestor;
        genelist[i].recentancestor=geneitems[i].recentancestor;*/
    }

    return nallspecies;
}
//.......................................................................................................................................................
void get_curgolgstats(uint64_t outgolgstats[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolgstats[ij] = golgstats[ij];                       // Note that golgstats is not dealt with in planes !
    }
}
//.......................................................................................................................................................
int get_sorted_popln_act( int gindices[], uint64_t genes[], int popln[], int activities[]) {
    int nspecies;
    int activitieshashx(int gindices[], uint64_t genes[], int popln[], int activities[]);
    nspecies=activitieshashx(gindices, genes, popln, activities);        // sets acttrace and returns current population arrays
    return(nspecies);
}
//-------------------------------------------------------------- comparison fns -------------------------------------------------------------------------
int cmpfunc (const void * pa, const void * pb) {
   // return ( *(int*)pa - *(int*)pb );
   return ((*(const uint64_t *)pa > *(const uint64_t *)pb)  ? 1 : -1);
}
//.......................................................................................................................................................
int cmpfunc1 ( const void *pa, const void *pb ) {
    const uint64_t *a = (const uint64_t *) pa;
    const uint64_t *b = (const uint64_t *) pb;
    if(a[1] == b[1])
        return a[0] > b[0] ? 1 : -1;
    else
        return (int) (b[1] - a[1]);
}
//.......................................................................................................................................................
int cmpfunc2 (const void * pa, const void * pb) {
    return ( genotypes[*(const int*)pa] > genotypes[*(const int*)pb] ? 1 : -1);
}
//.......................................................................................................................................................
int cmpfunc3 (const void * pa, const void * pb) {
    return ( geneitems[*(const int*)pa].popcount < geneitems[*(const int*)pb].popcount ? 1 : -1);
}
//.......................................................................................................................................................
int cmpfunc3c (const void * pa, const void * pb) {
    return ( cloneitems[*(const int*)pa].popln < cloneitems[*(const int*)pb].popln ? 1 : -1);
}
//.......................................................................................................................................................
int cmpfunc3q (const void * pa, const void * pb) {
    return ( quaditems[*(const int*)pa].pop1s < quaditems[*(const int*)pb].pop1s ? 1 : -1);
}
//.......................................................................................................................................................
int cmpfunc3qs (const void * pa, const void * pb) {
    uint64_t a,b;
    int na,nb;
    a = (uint64_t) *(const int*)pa;
    b = (uint64_t) *(const int*)pb;
    POPCOUNT64C(a,na);
    POPCOUNT64C(b,nb);
    return ( na < nb ? 1 : -1);
}
//.......................................................................................................................................................
int cmpfunc4 (const void * pa, const void * pb) {
   return ( geneitems[*(const int *)pa].firsttime > geneitems[*(const int *)pb].firsttime ? 1 : -1);
}
//.......................................................................................................................................................
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
//.......................................................................................................................................................
int cmpfunc5c (const void * pa, const void * pb) {               // sort according to ancestry in clonealogytrace

   int i1,i2,ij1,ij2,j;
   uint64_t birthid1,birthid2;
   i1=*(const int *)pa; i2=*(const int *)pb;

   for (j=0;j<clonealogydepth;j++) {
        ij1 = i1+j*N; ij2 = i2+j*N;
        birthid1=working[ij1]; birthid2=working[ij2];
        if(birthid1!=birthid2) return((((birthid1 > birthid2) && (birthid1!=rootclone)) || (birthid2==rootclone)) ? 1 : -1);
    }
    return(0);
}
//.......................................................................................................................................................
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
//.......................................................................................................................................................
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
//----------------------------------------------------------------- stats routines ----------------------------------------------------------------------
void countconfigs() {        // count translated configs specified by offset array
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
//.......................................................................................................................................................
void tracestats(uint64_t gol[],uint64_t golg[], uint64_t golgstats[], int NN2) { // trace various stats over time of the simulation
    int ij,cnt,k,d,dc,gt[4],st[10];
    uint64_t gene,statflag;

    if (!(diagnostics & diag_general_statistics)) {
        fprintf(stderr,"statistics collection not enabled in C\n");
        return;
    }
 
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
//-------------------------------------------------------------- countspecies ---------------------------------------------------------------------------
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
        else if ((selection == 2)||(selection == 3)) {                      // cyclic 4 species model
            fitness = nones&0x3;                                            // fitness is species class
        }

        // fprintf(stderr,"count species %d with gene %llx has counts %llu and %d ones, fitness %llu\n",k, golgsc[k][0],golgsc[k][1],nones,fitness);
    }
    //fprintf(stderr,"rulemod\trepscheme\tselection\toverwritemask\tsurvival\n");
    //fprintf(stderr,"%d\t%d\t\t%d\t\t%d\t\t%d\n",rulemod,repscheme,selection,overwritemask,survivalmask);
    //fprintf(stderr,"pmutmask\tinit1\tinitr\tncoding\tstartchoice\n");
    //fprintf(stderr,"%x\t\t%d\t%d\t%d\t%d\n",pmutmask,initial1density,initialrdensity,ncoding,startgenechoice);
}
//.......................................................................................................................................................
void countspecies() {                                                       // counts current species without using hash tables
    countspecies1(gol, golg, N2);
}
//.......................................................................................................................................................
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
//------------------------------------------------------------ totalpoptrace ---------------------------------------------------------------------------
int totalpoptrace(uint64_t gol[]) {   /* calculates and returns current population size and store in scrolling population trace array npopulation */
    int x, i, ij, npop;
    
    if (totdisp>=N) {                                               // 1 pixel to left scroll of population when full
            for(i=0;i<N;i++) npopulation[i]=npopulation[(i+1)&Nmask];
            x = N-1;
    }
    else x = totdisp;
    
    for(npop=ij=0;ij<N2;ij++) {                                     // calculate current population size
        npop+=(gol[ij]&0x1L) ? 1 : 0;
    }
    npopulation[x]=npop;
    
    return npop;
}

//------------------------------------------------------------ activitieshash ---------------------------------------------------------------------------
extern inline unsigned int genefnindex( uint64_t gene, uint64_t mask, int indexoff[]) {
    int k,d;
    uint64_t g,gf;
    uint64_t indexmask = 0xffffff;
    POPCOUNT64C(mask, d)
    if (d > 24) d=24;                                           // max of first 24 of active lut sites in mask allowed into genefnactivities array
    g = gene&mask;
    gf = 0ull;
    for (k=d;k>=0;k--) {
        gf |= (g>>indexoff[k])&0x1ull;
        gf <<= 1;
    }
    gf &= indexmask;
    return((unsigned int) gf);
}
//.......................................................................................................................................................
int activitieshash() {  /* count activities of all currently active gene species */
    int i, j, jmax, k, ij, ij1, x, nchist, nrhist, nspecies, nspeciesnow, cnt0, cnt1;
    int *gindices,*popln,*activities;
    double act,pop;
    uint64_t *genes, *traceptr;
    uint64_t gene,sbmask;
    const int maxact = 10000;
    int indexoff[24];
    
    sbmask = 0ull;
    if (activityfnlut) {
        for (j=0;j<24;j++) indexoff[j] = 0;
        sbmask = (((uint64_t) birthmask) << 32) | (uint64_t) survivalmask;
        for (j=k=0;k<64;k++) {
            if (j>=24) break;
            if ((sbmask>>k)&0x1ull) indexoff[j++]=k;
        }
        // fprintf(stderr,"indexoff:");for (j=0;j<24;j++) fprintf(stderr," %4d",indexoff[j]);fprintf(stderr,"\n");
    }

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );

    // if(gindices != NULL && col) free(gindices);
    for (i=0,nspeciesnow=0; i<nspecies; i++)
        nspeciesnow+=geneitems[i].popcount ? 1 : 0;

    gindices = (int *) malloc(nspeciesnow*sizeof(int));

    for (i=j=0; i<nspecies; i++) {
        if(geneitems[i].popcount) {
            gindices[j]=i;                                           // if col is 0 then the array gindices must be passed with sufficient length
            j++;
        }
        //gindices[j]=(geneitems[i].popcount) ? i : gindices[j--];   // if col is 0 then the array gindices must be passed with sufficient length
        //j++;
    }

    if(activityfnlut) {
        genes = (uint64_t *) malloc(nspeciesnow*sizeof(uint64_t));
        popln = (int *) malloc(nspeciesnow*sizeof(int));
        activities = (int *) malloc(nspeciesnow*sizeof(int));
        memset(popln,0,sizeof(int)*nspeciesnow);
        memset(activities,0,sizeof(int)*nspeciesnow);
        memset(genefnindices,0,sizeof(int)*(1<<24));
        genes = (uint64_t *) malloc(nspeciesnow*sizeof(uint64_t));
        for (i=0,jmax=1; i<nspeciesnow; i++) {
            k=genefnindex( genotypes[gindices[i]], sbmask, indexoff);
            j=genefnindices[k];
            if(!j) j=genefnindices[k]=jmax++;
            // genes[j-1]=genotypes[gindices[i]]&sbmask;
            genes[j-1]=genotypes[gindices[i]];
            popln[j-1]+=geneitems[gindices[i]].popcount;
            activities[j-1]+=geneitems[gindices[i]].activity;
        }
        nspeciesnow = jmax-1;
    }
    else {
        if (nspeciesnow > maxact) {                                      //sort in order of decreasing population
            qsort(gindices, nspeciesnow, sizeof(int), cmpfunc3);         // sort in decreasing count order
            nspeciesnow = maxact;
        }
        genes = (uint64_t *) malloc(nspeciesnow*sizeof(uint64_t));
        popln = (int *) malloc(nspeciesnow*sizeof(int));
        activities = (int *) malloc(nspeciesnow*sizeof(int));
        for (i=0; i<nspeciesnow; i++) {
            genes[i]=genotypes[gindices[i]];
            popln[i]=geneitems[gindices[i]].popcount;
            activities[i]=geneitems[gindices[i]].activity;
        }
    }

    if (totdisp>=N) {                                               // 1 pixel to left scroll when full
        for(ij=0;ij<N2;ij++) {
            ij1 = ((ij+1)&Nmask)+((ij>>log2N)<<log2N);              // (i+1)%N+j*N;
            // if(ij1>=N2) fprintf(stderr,"error in scroll of acttrace\n");
            acttrace[ij]=acttrace[ij1];
            poptrace[ij]=poptrace[ij1];
        }
        x=N-1;
    }
    else x=totdisp;

    for(i=0;i<N;i++) acttrace[x+i*N]=rootgene;                      // set column gray
    for(i=0;i<N;i++) poptrace[x+i*N]=rootgene;                      // set column gray
    //for(i=ymax1=0;i<nspeciesnow;i++)                              // ymax1 is current maximum of activities
    //    ymax1 = (activities[i]>ymax1) ? activities[i] : ymax1;
    // if (ymax1>ymax) ymax = ymax*2;                               // autoscale of activities
    // if (ymax1<ymax/2) ymax = ymax/2;                             // autoscale of activities
    for(j=0;j<nspeciesnow;j++) {
        gene = genes[j];
                                                                            // rescale populations and activities with saturation
        act = (double) activities[j];
        pop = (double) popln[j];
        // activities[j] = N-1 - (activities[j] * (N-1)) / ymax;    // linear scale, needs truncation if ymax superceded
        // activities[j] = (N-1) - (int) ((N-1)*log2(act)/log2ymax);// logarithmic scale, suffers from discrete steps at bottom
        activities[j] = (N-1) - (int) ((N-1)*act/(act+(double)ymax));
        popln[j] = (N-1) - (int) ((N-1)*pop/(pop+(double)ymax/10.));
        
        ij = (x&Nmask)+activities[j]*N;
        if(acttrace[ij]==rootgene)                                  // only one genotype to plot
            acttrace[ij] = gene;
        else {                                                      // plot species color with largest current population size if choice of multiple
            if(activityfnlut) {
                if((k=genefnindices[genefnindex( acttrace[ij], sbmask, indexoff)])) cnt0 = popln[k-1];
                else cnt0 = 0;
                cnt1 = popln[j];
                if(cnt1 >= cnt0) acttrace[ij]=gene;
            }
            else {
                if((genedataptr = (genedata *) hashtable_find(&genetable, acttrace[ij])) != NULL) cnt0 = genedataptr->popcount;
                else cnt0 = 0;
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) cnt1 = genedataptr->popcount;
                else cnt1 = 0;
                if(cnt1 >= cnt0) acttrace[ij]=gene;
            }
        }
        
        ij = (x&Nmask)+popln[j]*N;
        if(poptrace[ij]==rootgene)                                  // only one genotype to plot
            poptrace[ij] = gene;
        else {                                                      // plot species color with largest current activity if choice of multiple
            if(activityfnlut) {
                if((k=genefnindices[genefnindex( acttrace[ij], sbmask, indexoff)])) cnt0 = activities[k-1];
                else cnt0 = 0;
                cnt1 = activities[j];
                if(cnt1 >= cnt0) acttrace[ij]=gene;
            }
            else {
                if((genedataptr = (genedata *) hashtable_find(&genetable, poptrace[ij])) != NULL)
                   cnt0 = genedataptr->activity;
                else cnt0 = 0;
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL)
                   cnt1 = genedataptr->activity;
                else cnt1 = 0;
                if(cnt1 >= cnt0) poptrace[ij]=gene;
            }
        }
    }
    if(totdisp<N*nNhist) {
        nchist=totdisp/N; nrhist=totdisp-nchist*N;
        traceptr=&acttrace1[N2*nchist];
        for(j=0;j<N;j++) traceptr[nrhist+j*N]=acttrace[x+j*N];
        traceptr=&poptrace1[N2*nchist];
        for(j=0;j<N;j++) traceptr[nrhist+j*N]=poptrace[x+j*N];
        npopulation1[totdisp]=npopulation[x];
    }
    free(gindices);free(activities);free(genes);free(popln);
    return(nspeciesnow);
}
//.......................................................................................................................................................
int activitieshashx(int gindices[], uint64_t genes[], int popln[], int activities[]) {  /* python interface to count activities of all currently active species, no display */
    int i, j, nspecies, nspeciesnow;
    const int maxact = 10000;

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );

    // if(gindices != NULL && col) free(gindices);
    for (i=0,nspeciesnow=0; i<nspecies; i++)
        nspeciesnow+= geneitems[i].popcount ? 1 : 0;

    if (nspeciesnow>10000) return(-1);                              // exit with error need to allocate more space in python
    for (i=j=0; i<nspecies; i++) {
        gindices[j]=geneitems[i].popcount ? i : gindices[j--];      // if col is 0 then the array gindices must be passed with sufficient length
        j++;
    }
    qsort(gindices, nspeciesnow, sizeof(int), cmpfunc3);            // sort in decreasing count order

    if (nspeciesnow > maxact) nspeciesnow = maxact;

    for (i=0; i<nspeciesnow; i++) {
        genes[i]=genotypes[gindices[i]];
        popln[i]=geneitems[gindices[i]].popcount;
        activities[i]=geneitems[gindices[i]].activity;
    }

    return(nspeciesnow);                                            // exit here without doing display

}
//.......................................................................................................................................................
int activitieshashquad() {  /* count activities of all currently active quad images of connected components */
    int i, j, ij, ij1, x, nspecies, nspeciesnow, popcnt0, popcnt1;
    int *qindices,*popln,*activities;
    int qsindices[65536],qsallindices[65536];
    quadnode *q;
    double act;
    uint64_t *qids;                                                   // 64 bit ids for components used for colouring and ID
    uint64_t qid;
    const int maxact = 10000;

    nspecies = hashtable_count(&quadtable);
    quadkeys = hashtable_keys(&quadtable);
    quaditems = (quadnode*) hashtable_items( &quadtable );


    if(nhistG && (totsteps%nhistG == 0)) {                            // collect cumulative pattern size histograms of entire hash table
        for(i=0;i<=log2N;i++) histcumlogpattsize[i]=0;
        for(i=0;i<=N;i++) histcumpixelssqrt[i]=0;
        for (i=0; i<nspecies; i++) {
            histcumlogpattsize[log2upper(quaditems[i].size)]++;
            histcumpixelssqrt[sqrtupper(quaditems[i].pop1s)]++;
        }
    }
    
    if (totdisp>=N) {                                                 // 1 pixel to left scroll when full
        for(ij=0;ij<N2;ij++) {
            ij1 = ((ij+1)&Nmask)+((ij>>log2N)<<log2N);                // (i+1)%N+j*N;
            // if(ij1>=N2) fprintf(stderr,"error in scroll of acttraceq\n");
            acttraceq[ij]=acttraceq[ij1];
            acttraceqt[ij]=acttraceqt[ij1];
        }
        x=N-1;
    }
    else x=totdisp;
    for(i=0;i<N;i++) acttraceq[x+i*N]=rootgene;                       // set column gray, rootgene is used as unique pattern mapped to gray as for gene activities
    
    for (i=0,nspeciesnow=0; i<nspecies; i++)
        nspeciesnow+=quaditems[i].lasttime==totsteps ? 1 : 0;

    if (nspeciesnow) {
        qindices = (int *) malloc(nspeciesnow*sizeof(int));

        for (i=j=0; i<nspecies; i++) {
            if(quaditems[i].lasttime==totsteps) qindices[j++]=i;      // if col is 0 then the array qindices must be passed with sufficient length
        }
        if (nspeciesnow > maxact) {                                   // sort in order of decreasing pixel count
            qsort(qindices, nspeciesnow, sizeof(int), cmpfunc3q);
        }
        if (nspeciesnow > maxact) nspeciesnow = maxact;

        qids = (uint64_t *) malloc(nspeciesnow*sizeof(uint64_t));     // allocate arrays
        popln = (int *) malloc(nspeciesnow*sizeof(int));
        activities = (int *) malloc(nspeciesnow*sizeof(int));

        for (i=0; i<nspeciesnow; i++) {                               // set arrays of ids, popln (nr 1 pixels), and activities from hash table
            qids[i]=quadkeys[qindices[i]];
            popln[i]=quaditems[qindices[i]].pop1s;
            activities[i]=quaditems[qindices[i]].topactivity;
        }

        for(j=0;j<nspeciesnow;j++) {                                 // main loop to construct new display column for activities
            act = (double) activities[j];
            activities[j] = (N-1) - (int) ((N-1)*act/(act+(double)ymaxq));
            qid = qids[j];
            ij = (x&Nmask)+activities[j]*N;
            if(acttraceq[ij]==rootgene) {                            // first quadtype to plot at this position
                acttraceq[ij] = qid;
                acttraceqt[ij] = 1;                                  // tye of entry is quadtree (not small pattern)
            }
            else {                                                   // plot species color with largest current pop1s size if choice of multiple
                if((q = (quadnode *) hashtable_find(&quadtable, acttraceq[ij])) != NULL)
                   popcnt0 = q->pop1s;
                else popcnt0 = 0;
                if((q = (quadnode *) hashtable_find(&quadtable, qid)) != NULL)
                   popcnt1 = q->pop1s;
                else popcnt1 = 0;
                if(popcnt1 >= popcnt0) {
                    acttraceq[ij]=qid;
                    acttraceqt[ij] = 1;
                }
            }
        }
        free(qindices);free(activities);free(qids);free(popln);
    }
    
    if (nspeciesnow<maxact) {                                        // overlay activities of smallpatts up to maxact
        for (nallspeciessmall=nspeciessmall=i=0;i<65536; i++) {
            if (smallpatts[i].topactivity) {
                qsallindices[nallspeciessmall++]=i;
                if (smallpatts[i].lasttime == totsteps) qsindices[nspeciessmall++]=i;   // indices of current patterns
            }
        }
        if (nspeciesnow+nspeciessmall > maxact) {                    //sort in order of decreasing pixel count
            qsort(qsindices, nspeciessmall, sizeof(int), cmpfunc3qs);
            nspeciessmall = maxact-nspeciesnow;
        }
        
        activities = (int *) malloc(nspeciessmall*sizeof(int));
        for (i=0; i<nspeciessmall; i++) {
            activities[i]=smallpatts[qsindices[i]].topactivity;
        }
        
        for(j=0;j<nspeciessmall;j++) {
            act = (double) activities[j];
            activities[j] = (N-1) - (int) ((N-1)*act/(act+(double)ymaxq));
            qid = qsindices[j];
            ij = (x&Nmask)+activities[j]*N;
            if(acttraceq[ij]==rootgene) {                          // first quadtype to plot at this position (prefer quad over smallpatts)
                acttraceq[ij] = qid;
                acttraceqt[ij] = 0;
            }
            else if (!acttraceqt[ij]) {                            // plot species color with largest current pop1s size if choice of multiple
                POPCOUNT64C(((uint64_t) acttraceq[ij]), popcnt0);
                POPCOUNT64C(((uint64_t) qid), popcnt1);
                if(popcnt1 >= popcnt0) {
                    acttraceq[ij]=qid;
                    acttraceqt[ij] = 0;
                }
            }
        }
        free(activities);
        if(totdisp<N*nNhist) {
            for(i=0;i<N;i++) acttraceq1[totdisp+i*N]=acttraceq[x+i*N];
            for(i=0;i<N;i++) acttraceqt1[totdisp+i*N]=acttraceqt[x+i*N];
        }
    }
    
    return(nspeciesnow);
}
//--------------------------------------------------------------- genealogies ---------------------------------------------------------------------------
int get_genealogies(genedata genealogydat[], int narraysize) {  /* genealogies of all currently active species */
    int j, jmax, i, ij, k, nspecies, nspeciesnow;
    unsigned int birthstep;
    int j1, j2, j3, activity, gorder[N];
    uint64_t gene, ancgene, nextgene, genealogy1[N];
    uint64_t *curgen,*curgenealogy;
    int *gindices,*popln,*activities;
    genedata genedummy = {0,0,0,0,0,0,0,0ull,rootgene,rootgene};  // default data structure for gene data: note I put 4th member lastextinction -1 which is not a valid timestep

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );

    for (i=nspeciesnow=0; i<nspecies; i++)
        if(geneitems[i].popcount) nspeciesnow++;

    gindices = (int *) malloc(nspecies*sizeof(int));
    for (i=j=0; i<nspecies; i++) {
        if(geneitems[i].popcount) {
            gindices[j]=i;
            j++;
        }
        else gindices[nspeciesnow+i-j]=i;
    }

    if(!narraysize) {
        qsort(gindices, nspeciesnow, sizeof(int), cmpfunc3);// sort in decreasing population size order
        // qsort(gindices, nspeciesnow, sizeof(int), cmpfunc4);// alternatively, sort in increasing birthstep order
        if (nspeciesnow>N) nspeciesnow=N;                           // can only display at most N species, chose oldest
        curgenealogy=genealogy1;
    }
    else curgenealogy=working;

    popln = (int *) malloc(nspeciesnow*sizeof(int));
    activities = (int *) malloc(nspeciesnow*sizeof(int));

    for (i=0; i<nspeciesnow; i++) {
        popln[i]=geneitems[gindices[i]].popcount;
        activities[i]=geneitems[gindices[i]].activity;
    }

    if(narraysize) {                                                // need to allocate data for genealogy array in python
        for (i=jmax=0; i<nspeciesnow; i++) {                        // calculate max depth in genealogy jmax
            gene=genotypes[gindices[i]];
            ancgene=geneitems[gindices[i]].firstancestor;
            working[0]=gene;                                        // use working instead of genealogy1 array here as length needed may be larger than N
            for (j=1;;j++) {
                gene=ancgene;
                if(gene==rootgene) break;                           // reached root, exit j loop
                else {
                    for (k=0;k<j;k++) if (gene==working[k]) {gene=generepeat;break;};   // if gene already in ancestry, break with generepeat
                    if(gene==generepeat) { j=j+1;break;}
                    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
                        ancgene=genedataptr->firstancestor;
                    }
                    else fprintf(stderr,"ancestor not found in genealogies\n");
                }
                working[j]=gene;
            }
            if (j>jmax) jmax=j;
        }
        
        if(narraysize < (jmax+1)*nspeciesnow){
            fprintf(stderr,"get_genealogies(): narraysize not large enough.  Must be at least %d\n",nspeciesnow*(jmax+1));
            return(-1);
        }
        curgen = (uint64_t *) calloc(jmax,sizeof(uint64_t)); // current genealogy array
    }
    else {                                                        // display version
        for(ij=0;ij<N2;ij++) working[ij]=rootgene;                // set field to rootgene as background
        jmax=N-1;
        curgen = NULL;
    }
    
    activitymax=0;
    for (i=ij=genealogydepth=0; i<nspeciesnow; i++) {
        gene=genotypes[gindices[i]];                              // do not need to copy array to genes since only needed here
        if(narraysize) curgen[0]=gene;
        ancgene=geneitems[gindices[i]].firstancestor;
        activity=geneitems[gindices[i]].activity;
        if(activity>activitymax) activitymax=activity;
        curgenealogy[0]=gene;
        if(!narraysize) working[i]=gene;                          // ij = i for j=0
        for (j=k=1;j<=jmax;j++) {                                 // go back at most jmax links in genealogy
            gene=ancgene;
            if(gene==rootgene) break;                             // reached root, exit j loop
            else {
                for (k=0;k<j;k++) if (gene==curgenealogy[k]) {gene=generepeat;break;};  // if gene already in ancestry, break with generepeat
                if(gene==generepeat) { if(narraysize) curgen[j]=gene;else working[i+j*N]=gene;j=j+1;break;}
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
                    ancgene=genedataptr->firstancestor;
                    activity = genedataptr->activity;
                    if(activity>activitymax) activitymax=activity;
                }
                else fprintf(stderr,"ancestor not found in genealogies\n");
            }
            curgenealogy[j]=gene;
            if(narraysize) curgen[j]=gene; else working[i+j*N]=gene;
        }
        if (j>genealogydepth) genealogydepth=j;
        if(narraysize) {                                        // copy genealogy to python
            for(k=0; k<j; k++){
                if((genedataptr = (genedata *) hashtable_find(&genetable, curgen[k])) != NULL) {
                    genealogydat[ij+k] = *genedataptr; // cf comment on shallow copy in get_genes()
                }
                else fprintf(stderr,"ancestor not found in genealogies\n");
            }
            for(k=j; k<=jmax; k++){                             // fill list up to uniform jmax with dummy data
                genealogydat[ij+k] = genedummy;
            }
            ij += (jmax+1);
        }
    }
    if(narraysize) {                                            // if narraysize then return from routine here
        free(curgen);
        genealogydepth = jmax+1;
        return(genealogydepth);
    }
                                                                // reverse ancestries to allow comparison at same number of speciations
    for (i=0; i<nspeciesnow; i++) {
        for(j=0;j<N;j++) {
            if (working[i+j*N]==rootgene ) break;  // || working[i+j*N]==generepeat
        }
        for(j1=0;j1<(j>>1);j1++) {
            gene=working[i+(j-j1-1)*N];
            working[i+(j-j1-1)*N]=working[i+j1*N];
            working[i+j1*N]=gene;
        }
    }
    for (i=0; i<N; i++) gorder[i]=i;
    qsort(gorder, nspeciesnow, sizeof(int), cmpfunc5);          // sort according to ancestral lines - use cmpfunc5 to sorting genes laterally via gene value
    //qsort(gorder, nspeciesnow, sizeof(int), cmpfunc6);        // sort according to ancestral lines - use cmpfunc6 to sorting genes laterally via activity
    //qsort(gorder, nspeciesnow, sizeof(int), cmpfunc7);        // sort according to ancestral lines - use cmpfunc7 to sort genes laterally via population size

    for (i=0;i<N;i++) if((gorder[i]<0)||(gorder[i]>=N)) fprintf(stderr,"step %d error in gorder out of bounds at i = %d with value %d\n",totsteps,i,gorder[i]);

    for(ij=0;ij<N2;ij++) genealogytrace[ij]=rootgene;           // initialize genealogytrace to root gene before drawing part of it

    if(colorfunction==7 || colorfunction2==7) {                 // time trace of genealogies: NB simultaneus display of 6 and 7 not possible :> 2x 7 or 2x 6
      birthstep=0;
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
                    if((genedataptr = (genedata *) hashtable_find(&genetable, nextgene)) != NULL) {
                        birthstep = (unsigned int) genedataptr->firsttime;
                    }
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
    else {                                                      // species changes only trace (colorfunction == 6)
      for(i=0;i<nspeciesnow;i++) {
        for(j=0;j<jmax;j++) {
            ij=i+j*N;
            genealogytrace[ij]=working[gorder[i]+j*N];
        }
      }
    }
    free(gindices);free(activities);free(popln);
    genealogydepth = jmax+1;

    return(genealogydepth);
}
//--------------------------------------------------------------- clonealogies --------------------------------------------------------------------------
int clonealogies() {                                            // genealogies of all clones
    int j, jmax, i, ij, nclones, nclonesnow, birthstep;
    int j1, j2, j3, activity, gorder[N];
    uint64_t birthid, parentid, nextclone;
    // uint64_t gene, ancgene, nextgene;
    int *gindices,*popln,*activities,*birthsteps;

    nclones = hashtable_count(&clonetable);
    clones = hashtable_keys(&clonetable);
    cloneitems = (clonedata*) hashtable_items( &clonetable );

    for (i=nclonesnow=0; i<nclones; i++)
        if(cloneitems[i].popln) nclonesnow++;

    gindices = (int *) malloc(nclones*sizeof(int));
    for (i=j=0; i<nclones; i++) {
        if(cloneitems[i].popln) {
            gindices[j]=i;
            j++;
        }
        else gindices[nclonesnow+i-j]=i;
    }

    qsort(gindices, nclonesnow, sizeof(int), cmpfunc3c);        // sort in decreasing population size order

    if (nclonesnow>N) nclonesnow=N;                             // can only display at most N clones, chose oldest

    popln = (int *) malloc(nclonesnow*sizeof(int));
    activities = (int *) malloc(nclonesnow*sizeof(int));
    birthsteps = (int *) malloc(nclonesnow*sizeof(int));

    for (i=0; i<nclonesnow; i++) {
        popln[i]=cloneitems[gindices[i]].popln;
        activities[i]=cloneitems[gindices[i]].activity;
        birthsteps[i]=(unsigned int) ((cloneitems[gindices[i]].birthid)>>32);
    }

    for(ij=0;ij<N2;ij++) working[ij]=rootclone;                 // set field to rootclone as background
    activitymax=0;
    for (i=jmax=0; i<nclonesnow; i++) {
        birthid=clones[gindices[i]];                            // do not need to copy array since only needed here
        parentid=cloneitems[gindices[i]].parentid;
        activity=cloneitems[gindices[i]].activity;
        if(activity>activitymax) activitymax=activity;
        working[i]=birthid;                                     // ij = i for j=0
        for (j=1;j<N;j++) {                                     // go back at most N links in clonealogy
            birthid=parentid;
            if(birthid&rootclone) break;                       // reached root, exit j loop
            else {
                if((clonedataptr = (clonedata *) hashtable_find(&clonetable, birthid)) != NULL) {
                    parentid=clonedataptr->parentid;
                    activity = clonedataptr->activity;
                    if(activity>activitymax) activitymax=activity;
                }
                else fprintf(stderr,"ancestor not found in clonealogies\n");
            }
            ij = i+j*N;
            working[ij]=birthid;
        }
        if (j>jmax) jmax=j;
    }
    clonealogydepth = jmax;

                                                                //reverse ancestries to allow comparison at same number of clonal speciations
    for (i=0; i<nclonesnow; i++) {
        for(j=0;j<N;j++) {
            if (working[i+j*N]&rootclone) break;
        }
        for(j1=0;j1<(j>>1);j1++) {
            birthid=working[i+(j-j1-1)*N];
            working[i+(j-j1-1)*N]=working[i+j1*N];
            working[i+j1*N]=birthid;
        }
    }
    for (i=0; i<N; i++) gorder[i]=i;
    qsort(gorder, nclonesnow, sizeof(int), cmpfunc5c);          // sort according to ancestral lines - use cmpfunc5 to sorting clones laterally via clone value

    for (i=0;i<N;i++) if((gorder[i]<0)||(gorder[i]>=N)) fprintf(stderr,"step %d error in gorder out of bounds at i = %d with value %d\n",totsteps,i,gorder[i]);

    for(ij=0;ij<N2;ij++) clonealogytrace[ij]=rootclone;          // initialize clonealogytrace to root clone before drawing part of it

    if(colorfunction==7 || colorfunction2==7) {                 // time trace of clonealogies (takes precedence with two difft displays
      birthstep=0;
      for(i=0;i<nclonesnow;i++) {
        for(j=0,j1=0;j<jmax;j++) {
            if(gorder[i]>=nclonesnow) fprintf(stderr,"error in clonealogies gorder at i=%d, order value %d out of range\n",i,gorder[i]);
            ij = gorder[i]+j*N;
            birthid = working[ij];
            ij+=N;
            if(ij<N2) {
                nextclone = working[ij];
                if(nextclone&rootclone) birthstep=totsteps;
                else {
                    if((clonedataptr = (clonedata *) hashtable_find(&clonetable, nextclone)) != NULL) birthstep = (unsigned int) (clonedataptr->birthid>>32);
                    else fprintf(stderr,"ancestor %llx not found at (%d,%d) in clonealogies during birthstep extraction\n",nextclone,ij&Nmask,ij>>log2N);
                }
            }
            else birthstep=totsteps;
            j2 = birthstep*N/totsteps;
            for (j3=j1;j3<j2;j3++) {
                ij = i+j3*N;
                clonealogytrace[ij]=birthid;
            }
            j1 = j2;
        }
      }
      // for(i=nclonesnow;i<N;i++) for(j=0;j<N;j++) clonealogytrace[gorder[i]+j*N]=rootclone;
    }
    else {                                                      // species changes only trace (colorfunction == 6)
      for(i=0;i<nclonesnow;i++) {
        for(j=0;j<jmax;j++) {
            ij=i+j*N;
            clonealogytrace[ij]=working[gorder[i]+j*N];
        }
      }
    }
    free(gindices);free(activities);free(popln);free(birthsteps);
    return(jmax);
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------

void get_gliderinfo(uint64_t outgliderinfo[], int narraysize){               // put 7x7 pattern averaged match counts into outgliderinfo array
    uint64_t *gitmp, gene;
    int ij,k,nbhood,sum=0;
    unsigned int d1;
    if(narraysize!=51*8 && narraysize!=27*8 && narraysize!=11*8){
        fprintf(stderr,"get_gliderinfo():  wrong data size (should be array of 8*11,8*27 or 8*51)\n");
    }
    nbhood=narraysize>>3;                                     // divide by eight
    if (nbhood != it_nbhood*it_nbhood+2) {
        fprintf(stderr,"error: mismatch between C value of it_nbhood %d and python expectation %d\n",it_nbhood,nbhood);
        return;
    }
    gitmp=outgliderinfo;
    for (ij=0; ij<N2; ij++) {
        gene = golmix[ij];
        for (k=0;k<8;k++) {                                   // each direction: N E S W NE SE SW NW
            gitmp = outgliderinfo + k*nbhood;
            d1 = (int) ((gene>>(k<<3))&0xffull);              // differences for this direction
            if(d1==0xff)
                gitmp[nbhood-1]++;
            else{
                if(d1<0 || d1>nbhood-2){
                    fprintf(stderr, "get_gliderinfo:  bad difference count value:  %d.", d1);
                    return;
                }
                d1 = nbhood-2-d1;                             // change to matches
                gitmp[d1]++;
            }
        }
    }
    for(d1=0;d1<narraysize;d1++) sum += gitmp[d1];
    // fprintf(stderr,"in get_gliderinfo sum is %d\n",sum);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------
