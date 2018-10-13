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

const int log2N = 9;                // toroidal array of side length N = 2 to the power of log2N
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
int selection = 1;                  // fitness of 2 live neighbours: integer value (0), number of ones (1), Norman compete (2), 2 target coding (3)
int repscheme = 1;                  // replication scheme: 0 random choice , 1 XOR of all 3, 2 consensus, 3 unique det choice, 4 most different
int ncoding = 16;                   // maximal distances between sequences for fitness2 == 2
                                    // number of bits used to encode non zero bits in 2nd gene in sequence space for fitness2 == 3
int survival = 0x2;                 // survive mask for two (bit 1) and three (bit 0) live neighbours
int overwritemask = 0x2;            // bit mask for 4 cases of overwrite: bit 0. s==3  bit 1. special birth s==2
int initial1density = (1<<15)>>1;   // initial density of ones in gol as integer value, divide by 2^15 for true density
int initialrdensity = (1<<15)>>1;   // initial density of random genes in live sites, divide by 2^15 for true density
int startgenechoice = 8;            // selection for defined starting genes 0-8 (8 is random 0-7) otherwise choose only particular startgene
long unsigned int codingmask;       // mask for ncoding bits
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

    int k, kmin, nmut;
    int nb[8], ij, i, j, jp1, jm1, ip1, im1;
    long unsigned int s, gs, nb1i, randnr, randnr2, r2;
    long unsigned int nbmask, nbmaskr, nbmaskrm;
    long unsigned int newgene, livegenes[3],gdiff,gdiff0,gdiff1;
    long unsigned int s2or3, birth;
    int d0,d1,d2,d3,dd;
    long unsigned int gene2centre;
    int g0011,g0110,prey;
    
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
       s2or3 = (s>>2) ? 0L : (s>>1);                                       // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
       if (s2or3) {                                                        // if 2 or 3 neighbours alive
            if ((s<2)||(s>3)) fprintf(stderr,"s2or3 error s == %lu\n",s);
            birth = 0L;
            newgene = 0L;
            if (s&0x1L) {  // s==3                                                 // birth (with possible overwrite)
              if ((0x1L&overwritemask)|(0x1L&~gol[ij]) ) {                         // central site empty or overwrite mode
                birth = 1L;                                                        // birth flag
                if (repscheme & 0x1) {
                    for(k=0;k<s;k++) livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];  // live neighbour genes
                    if((livegenes[0]^livegenes[1])|(livegenes[0]^livegenes[2])) {  // genes not all same
                        gdiff  = livegenes[0]^livegenes[2];
                        gdiff0 = livegenes[0]^livegenes[1];
                        gdiff1 = livegenes[1]^livegenes[2];
                        if(gdiff0&&gdiff1&&gdiff) {                                //all three different
                            POPCOUNT64C(gdiff0,d0);
                            POPCOUNT64C(gdiff1,d1);
                            POPCOUNT64C(gdiff,d2);
                            d3 = d0!=d1? (d0!=d2 ? 4 : 1) : (d0!= d2 ? 2 : (d1!=d2 ? 3 : 0));
                        }
                        else if (gdiff1) livegenes[1]=livegenes[2];
                    }
                    else newgene= golg[nb[nb1i&0x7]];
                }
                else {
                  for (k=7,nbmask=0L;k>=0;k--) nbmask = (nbmask << 1) + gol[nb[k]];  // compute 8-bit mask of GoL states of 8 nbs, clockwise from top left
                  for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {                 // compute canonical rotation (minimum) of this mask
                    nbmaskr = ((nbmaskr & 0x1L)<<7) + (nbmaskr>>1);                // 8 bit rotate right
                    if (nbmaskr < nbmaskrm) {                                      // choose minimal value of mask rotation
                        nbmaskrm = nbmaskr;                                        // neighbor mask rotate min is current rotation
                        kmin = k;                                                  // no of times rotated to right
                    }
                  }

                  if ((repscheme>>1) & 0x1) newgene = golg[nb[kmin]];                // 2,3,6,7... deterministic choice of ancestor: replication of live neigbour in bit 0 of canonical pos
                  else {                                                             // repscheme = 0 or 4 for example
                    switch (nbmaskrm) {                //             x07   x0b   x13   x19   x0d   x15   x25
                        case 0x07 : k = 1; break;      // 00000111    012  <-
                        case 0x0b : k = 0; break;      // 00001011    ...   01.  <-
                        case 0x13 : k = 1; break;      // 00010011    ...   ..3   01.  <-
                        case 0x19 : k = 0; break;      // 00011001          ...   ...   0..  <-
                        case 0x0d : k = 3; break;      // 00001101                ..4   ..3   0.2  <-
                        case 0x15 : k = 2; break;      // 00010101                      ..4   ..3   0.2  <-
                        case 0x25 : k = 5; break;      // 00100101                            ...   ...   0.2  <-
                        default  : {                   //                                           ..4   ...
                                                       //                                                 .5.
                            fprintf(stderr,"Error in canonical rotation for three live neighbours \nnbmaskrm = %lx\n",nbmaskrm); k = 0;
                            fprintf(stderr,"Raw Neighbor Pattern: %lx No neighbors %lx\n",
                                nbmask, gol[nb[0]]+gol[nb[1]]+gol[nb[2]]+gol[nb[3]]+gol[nb[4]]+gol[nb[5]]+gol[nb[6]]+gol[nb[7]]);
                            fprintf(stderr,"\n");
                        } //case
                    } //switch
                    newgene = golg[nb[(kmin+k)&0x7]];                               // rotate unique nb k left (kmin) back to orig nb pat
                  }
                  //if (newgene == 0L) {
                  //    fprintf(stderr,"step %d Error with new gene zero: nbmask %lu nbmaskrm %lu kmin %d gol %lu golg %lx newgene %lx ij %d\n",totsteps,nbmask,nbmaskrm,kmin,gol[nb[kmin]],golg[nb[kmin]],newgene,ij);
                  //}
                }
              }
            }  // end else if s==3
            else {  // s==2                                                 // possible birth as exception to GoL rule
                if (rulemod) {                                              // special rule allowed if rulemod==1, NB no birth if all sequences same
                    if ((0x1L&(overwritemask>>1))|(0x1L&~gol[ij])) {        // either overwrite on for s==2 or central site is empty
                        for (k=0;k<s;k++)                                   // loop only over live neigbours
                            livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];    // live gene at neighbour site
                        if (selection==0) {                                 // use integer value of sequence as fitness
                            birth = (livegenes[0]^livegenes[1]) ? 1L: 0L;   // birth condition is two genes different
                            newgene = livegenes[0]>livegenes[1] ?  livegenes[0] : livegenes[1]; // choose one with larger gene to replicate
                        }
                        else {
                            POPCOUNT64C(livegenes[0],d0);
                            POPCOUNT64C(livegenes[1],d1);
                            if (selection==1) {                             // use number of ones in sequence as fitness
                                birth = (d0^d1) ? 1L: 0L;                   // birth condition if two genes different in number of ones
                                newgene= (d0>d1) ? livegenes[0] : livegenes[1];
                            }
                            else if (selection==2) {                        // use scissors-paper-stone-well game on number ones mod 4
                                d2=d0&0x3;
                                d3=d1&0x3;
                                birth = (d2^d3) ? 1L: 0L;                   // birth if 2 genes differ mod 4 in number of ones
                                newgene= 3-d2 ? livegenes[1] : (d2<d3 ? livegenes[1] : livegenes[0]);
                            }
                             else if (selection==3) {                       // birth if 2 genes differently functional (Not Yet Working)
                                gene2centre = (1L<<ncoding)-1L;             // first ncoding 1s in this sequence
                                gdiff  = livegenes[0]^livegenes[1];
                                gdiff0 = livegenes[0]^gene2centre;
                                gdiff1 = livegenes[1]^gene2centre;
                                POPCOUNT64C(gdiff,dd);
                                POPCOUNT64C(gdiff0,d2);
                                POPCOUNT64C(gdiff1,d3);
                                g0011 = d0<dd && d3<dd;
                                g0110 = d2<dd && d1<dd;
                                birth = (g0011 || g0110)  ? 1L: 0L;         // birth if 2 genes closer to two different targets than each other
                                newgene= g0011 ? ((d0<d3) ? livegenes[0] : livegenes[1]) : ((d2<d1) ? livegenes[0] : livegenes[1]);
                            }
                            else if (selection==4) {                        // birth if 2 genes obey 3 distance constraints < ncoding (NYW)
                                gdiff0 = ~livegenes[0];                     // find distance to 0x0 for gene 0
                                POPCOUNT64C(gdiff0,d0);                     // in d0
                                gdiff=livegenes[0]^livegenes[1];
                                POPCOUNT64C(gdiff,dd);
                                birth = (dd<ncoding && d0<ncoding && d1<ncoding) ? 1L: 0L; // birth if 2 genes close enough (< has higher priority than &&)
                                newgene= (d0>d1) ? livegenes[0] : (d0<d1 ? livegenes[1] : ((livegenes[0]>livegenes[1]) ?  livegenes[0] : livegenes[1]));
                            }
                            else if (selection==5) {                        // predator prey model : prey evolves to all 0, predator to all 1
                                prey = d0<32 || d1<32;                      // >=1 prey required for birth
                                gdiff=livegenes[0]^livegenes[1];
                                birth = (gdiff && prey) ? 1L: 0L;         // birth if different and >=1 prey)
                                prey = d0<32 && d1<32;                      // 2 prey : newgene is one with less ones, 1 prey : predator wins
                                newgene= (prey ? (d0<d1 ? livegenes[0] : livegenes[1]) : (d0<d1 ? livegenes[1] : livegenes[0]));
                            }
                            else fprintf(stderr,"Error: two live gene fitness value %d is not allowed\n",selection);
                        }
                    }
                }
            }

            if(birth){
                // if (gol[ij]) fprintf(stderr,"birth overwrite event ij %d newgene %lu s %lu\n",ij,newgene,s);
                RAND128P(randnr);                                               // inline exp so compiler recognizes auto-vec,
                // compute random events for single bit mutation, as well as mutation position nmut
                randnr2 = (randnr >> 24) & pmutmask;                                // extract bits from randnr for random trial for 0 on pmutmask
                r2 = randnr2?0L:1L;                                                   // 1 if lowest nlog2pmut bits of (bits 24-47 of randnr) are zero, else zero
                nmut = (randnr >> 48) & 0x3f;                                       // choose mutation position for length 64 gene : from bits 48:53 of randnr
                // complete calculation of newgol and newgolg, including mutation
                newgene = newgene ^ (r2<<nmut);                              // introduce single mutation with probability pmut = probmut
                newgol[ij]  =  1L;                                      // new game of life cell value: stays same or set to one from zero if birth
                newgolg[ij] =  newgene;                                      // if birth then newgene
            }
            else {
                if ((survival&s&0x1L)|((survival>>1)&(~s)&0x1L)|((~rulemod)&0x1L)) { // survival bit 0 and s==3, or (survival bit 1 and s==2) or not rulemod
                    newgol[ij]  = gol[ij];                                  // new game of life cell value same as old
                    newgolg[ij] = golg[ij];                                 // gene stays as before, live or not
                }
                else {
                    newgol[ij]  = 0L;                                       // new game of life cell value dead
                    newgolg[ij] = 0L;                                       // gene dies or stays dead
                }
            }
        }  // end if s2or3
        else {                                                              // else not birth or survival, 0 values for gol and gene
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

char *readFile(char *fileName)
{
  FILE *file;
  char *code = malloc(32* 32 * sizeof(char));
  file = fopen(fileName, "r");
  do
  {
    *code++ = (char)fgetc(file);

  } while(*code != EOF);
  fclose(file);
  return code;
}

int writeFile(char *fileName)
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
    int ij,ij1,k;
    long unsigned int g;
    long unsigned int *gol;
    long unsigned int *golg;
    static unsigned int rmask = (1 << 15) - 1;
    // Range: rand returns numbers in the range of [0, RAND_MAX ), and RAND_MAX is specified with a minimum value of 32,767. i.e. 15 bit

    long unsigned int startgenes[8];
    
    static int Nf = 0;
    char *golgin;
    
    srand(1234567);
    state[0] = rand();state[1] = rand();

    // writeFile("genepat.dat");

    rulemod = runparams[0];
    repscheme = runparams[1];
    selection = runparams[2];
    overwritemask = runparams[3];
    survival = runparams[4];

    nlog2pmut = simparams[0];
    pmutmask = (0x1 << nlog2pmut) - 1;
    initial1density = simparams[1];
    initialrdensity = simparams[2];
    ncoding = simparams[3];
    codingmask = (0x1L<<ncoding)-1;
    startgenechoice = simparams[4];
    
    fprintf(stdout,"runparams %d %d %d %d %d\n",runparams[0],runparams[1],runparams[2],runparams[3],runparams[4]);
    fprintf(stdout,"simparams %d %d %d %d %d\n",simparams[0],simparams[1],simparams[2],simparams[3],simparams[4]);

    startgenes[0] = 0xaaaa000000000000;
    startgenes[1] = 0xaaaa00000000ffff;
    startgenes[2] = 0xaaaa0000ffff0000;
    startgenes[3] = 0xaaaa0000ffffffff;
    startgenes[4] = 0xaaaaffff00000000;         // should have same functionality as 0
    startgenes[5] = 0xaaaaffff0000ffff;         //   ... as 1
    startgenes[6] = 0xaaaaffffffff0000;         //   ... as 2
    startgenes[7] = 0xaaaaffffffffffff;         //   ... as 3
    
    startgenes[0] = 0xffffffffffc00000;
    startgenes[1] = 0xffffffffffc00000;
    startgenes[2] = 0xffffffffffc00000;
    startgenes[3] = 0xffffffffffc00000;
    startgenes[4] = 0x00000fffffc00000;
    startgenes[5] = 0x00000fffffc00000;
    startgenes[6] = 0x00000fffffc00000;
    startgenes[7] = 0x00000fffffc00000;
  
    gol = planes[curPlane];
    golg = planesg[curPlane];
 
    if (Nf) {           // input from file
        golgin=readFile("genepat.dat");
        for (ij=0; ij<N2; ij++) {
            gol[ij] = 0;
            golg[ij] = 0;
        }
        for (ij1=0; ij1<32*32; ij1++) {
            ij=(N>>1)-16+(ij1&0x1f)+ N*((N>>1)-16+(ij1>>5));
            if (golgin[ij1] > 0)    {                   // if live cell
                gol[ij] = 1L;
                if(golgin[ij1] <= 8 ) golg[ij] = startgenes[golgin[ij1]-1];
                else golg[ij] = startgenes[7];
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
            }
            golg[ij] = g;
            // if (golg[ij] == 0L && gol[ij] != 0L) fprintf(stderr,"zero gene at %d\n",ij);
        }
        // for (ij=0; ij<40; ij++) fprintf(stderr,"gene at %d %lx\n",ij,golg[ij]);   // test first 40
    }
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

void countspecies(long unsigned int gol[], long unsigned int golg[], int params[], int N2, int nparams) {  /* counts numbers of all different species using qsort first */
    int ij, k, ijlast, nspecies, counts[N2], nones;
    long unsigned int last, golgs[N2], fitness;
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
    counts[k++]=N2-ijlast;
    nspecies = k;  // print including 0
    fprintf(stdout,"The number of different non-zero species is %d\n",nspecies);

    for (k=0,ij=0;k<nspecies;k++) {     // now condense array to give only different genes with counts
        // printf("species %4d with gene %x has counts %d\n",k, golgs[ij],counts[k]);
        golgs[k]=golgs[ij];
        ij = ij + counts[k];
    }

    for (k=0; k<nspecies; k++) { golgsc[k][0] = golgs[k];  golgsc[k][1] = counts[k];}  // initialize joint gene & count array
    qsort(golgsc, nspecies, sizeof(golgsc[0]), cmpfunc2);                   // sort in decreasing count order
    for (k=0; k<nspecies; k++) {
        last = golgsc[k][0];
        POPCOUNT64C(last, nones);
        fitness = 999;
        if (selection == 0) {                                               // 2-live neighbor fitness is integer value
	        fitness = last;
        }
        else if (selection == 1) {                                          // 2-live neighbor fitness is number of ones
	        fitness = (long unsigned) nones;
        }
        else if (selection == 2){                                           // non-neutral model based on presence of replicase gene
	        fitness = 999;                                                  // undefined, depends on competing sequence
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
        fprintf(stdout,"count species %d with gene %lx has counts %lu and %d ones, fitness %lu\n",k, golgsc[k][0],golgsc[k][1],nones,fitness);
    }
    fprintf(stdout,"cumulative activity = %lu\n",(N2 * (long unsigned int) totsteps) - emptysites);
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

int colorFunction = 0;

void colorgenes(long unsigned int gol[],long unsigned int golg[], int cgolg[], int N2) {
    long unsigned int gene, mask;
    int ij,d;
    int cnt;
    cnt=0;
    if(colorFunction){
	    for (ij=0; ij<N2; ij++) {
	        if (gol[ij]) {
		        gene = golg[ij];
                POPCOUNT64C(gene,d);
		        cgolg[ij] = 2 + 3*d;
	        }
	        else cgolg[ij] = 0;
	    }
    }
    else{
	    for (ij=0; ij<N2; ij++) {
            if (gol[ij]) {
                gene = golg[ij];
                if (gene == 0L) gene = 11778L; // random color for gene==0
                // mask = (gene * 11400714819323198549ul) >> (64 - 8);   // hash with optimal prime multiplicator down to 8 bits
                mask = (gene * 11400714819323198549ul) >> (64 - 32);   // hash with optimal prime multiplicator down to 32 bits
                cgolg[ij] = 1 + (int) mask;
                if (((int) mask) == 0) cnt++;
            }
            else cgolg[ij] = 0;
	    }
    }
    // fprintf(stderr,"hash count of zero genes %d\n",cnt);
}
