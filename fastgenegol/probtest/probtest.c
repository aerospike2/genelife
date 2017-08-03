//
//  maindispatch.c
//  probtest
//
//  Created by John McCaskill on 31.07.17.
//  Copyright © 2017 Ruhr Universität Bochum. All rights reserved.
//  This version uses the dispatch library to implement a concurrent loop
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

const int log2N = 7;                // toroidal array of side length N = 2 to the power of log2N
const int N = 0x1 << log2N;
const int N2 = N*N;                 // number of sites in toroidal array
const int Nmask = N - 1;            // bit mask for side length, used instead of modulo operation

                                    // Wikipedia "Xorshift" rewritten here as inline macro &
                                    // Vigna, Sebastiano. "xorshift*/xorshift+ generators and the PRNG shootout". Retrieved 2014-10-25.
static long unsigned int state[2]; // State for xorshift pseudorandom number generation. The state must be seeded so that it is not zero
inline long unsigned int rand128p() {
    long unsigned int x = state[0]; long unsigned int const y = state[1];
	state[0] = y;	x ^= x << 23;  state[1] = x ^ y ^ (x >> 17) ^ (y >> 26);
	return state[1] + y;}
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
inline void popcount2x32(long unsigned x, int *cnt1p, int *cnt2p) {
    x -= (x >> 1) & m1;             /* put count of each 2 bits into those 2 bits */
    x = (x & m2) + ((x >> 2) & m2); /* put count of each 4 bits into those 4 bits */
    x = (x + (x >> 4)) & m4;        /* put count of each 8 bits into those 8 bits */
    *cnt1p = ((x&hlo) * h02) >> 24;   /* left 8 bits of x + (x<<8) + (x<<16) + (x<<24) for x = x&hlo */
    *cnt2p = ((x>>32) * h02) >> 24;}  /* left 8 bits of y + (y<<8) + (y<<16) + (y<<24) for y = x>>32  */
#define POPCOUNT2X32(x, val1, val2) {      /* Separate Hamming weights of lower and upper 32 bits */  \
    x -= (x >> 1) & m1;             /* put count of each 2 bits into those 2 bits */ \
    x = (x & m2) + ((x >> 2) & m2); /* put count of each 4 bits into those 4 bits */ \
    x = (x + (x >> 4)) & m4;        /* put count of each 8 bits into those 8 bits */ \
    val1 = ((x&hlo) * h02) >> 24;   /* left 8 bits of x + (x<<8) + (x<<16) + (x<<24) for x = x&hlo */ \
    val2 = ((x>>32) * h02) >> 24;}  /* left 8 bits of y + (y<<8) + (y<<16) + (y<<24) for y = x>>32  */

const long unsigned int m3  = 0x0707070707070707; //for cumcount64c selects counter relevant bits only
#define CUMCOUNT64C(x, val) {       /* Assumes gene specifies 8 8-bit counters each with max value 7 */  \
    val = ((x & m3) * h01) >> 56;}  /* left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... */

void genelife_update (long unsigned int gol[], long unsigned int golg[], int ndsteps) {
	/* update GoL for toroidal field which has side length which is a binary power of 2 */
	/* encode without if structures for optimal vector treatment */
	int t, nones2, nones3, ij;
    // int i,j;
    // int ni2, ni3;
    // double xp2, xp3;
    long unsigned int randnr, r2, r3;
    unsigned int p2, p3;
    long unsigned int genef1;
    static double a2=1., a3=1.;
    static long unsigned int  newgol[N2],newgolg[N2];
    const double max32uint = 4294967295.;

    static int first = 1;
    
    if (first) {
        state[0] = rand();state[1] = rand();
        first = 0;
    }
    for (t=0; t<ndsteps; t++) {
	  for (ij=0; ij<N2; ij++) {                                               // loop over all sites of 2D torus with side length N
		//i = ij & Nmask;  j = ij >> log2N;                                   // row & column
        // randnr = rand128p();                                                   // expansion inline so compiler recognizes auto-vectorization options
        RAND128P(randnr);                                                     // expansion inline so compiler recognizes auto-vectorization options
        genef1 = golg[ij];                                                    // gene of central site
        //popcount2x32(genef1, &nones2, &nones3);
        POPCOUNT2X32(genef1, nones2, nones3);
        
        p2=(unsigned int) (max32uint*exp(-a2*nones2));                        // 15% 2s longer on 100000 than piecewise linear version below
        p3=(unsigned int) (max32uint*exp(-a3*nones3));
        //ni2= 1 + (int) (xp2=a2*nones2);
        //ni3= 1 + (int) (xp3=a3*nones3);
        //p2 = ((unsigned int)(max32uint*(1.-xp2+ni2))) >> ni2;               // piecewise linear reasonable approximation to exponential
        //p3 = ((unsigned int)(max32uint*(1.-xp3+ni3))) >> ni3;
        
        r2 = (randnr & hlo) < p2 ? 0x1L: 0;
        r3 = ((randnr>>32) & hlo) < p3 ? 0x1L: 0;
        newgol[ij] = r2+(r3<<1);
        newgolg[ij] = randnr;
      }
      
      for (ij=0; ij<N2; ij++) {
        gol[ij] = newgol[ij];        // copy new gol config to old one
        golg[ij] = newgolg[ij];      // copy new genes to old genes
      }
    }
}

void initialize_genes (long unsigned int golg[], long unsigned int gol[], int N2) {
	int ij,k;
    long unsigned int g;
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

int main(int argc, const char * argv[]) {
    int ij;
    static long unsigned int  gol[N2],golg[N2];
    printf("Starting...\n");
    for (ij=0; ij<N2; ij++) {
		gol[ij] = (rand() & 0x1)?1:0;
	}
    initialize_genes (golg, gol, N2);
    genelife_update (gol, golg, 100000);
    printf("Finished!\n");
    return 0;
}
