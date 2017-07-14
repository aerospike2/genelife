//
//  fastgol.c
//  fastggenegol
//
//  Created by John McCaskill on 14.07.17.
//  Copyright © 2017 European Center for Living Technology. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>

#define log2N 7

const int N = 0x1 << log2N;
const int N2 = N*N;
const unsigned int Nmask = N - 1;



void update (int gol[]) {
	/* update GoL for toroidal field which has side length which is a binary power of 2 */
	/* encode without if structures for optimal vector treatment */

	unsigned int ij, i, j, jN, jp1, jm1, ip1, im1, s;
	static unsigned int newgol[N2];

	for (ij=0; ij<N2; ij++) {
		i = ij & Nmask;
		j = ij >> log2N;
		jN = j*N;
		jp1 = ((j+1) & Nmask)*N;
		jm1 = ((j-1) & Nmask)*N;
		ip1 = (i+1) & Nmask;
		im1 = (i-1) & Nmask;
		s = gol[jN+ip1]+gol[jN+im1]+gol[jp1+i]+gol[jp1+ip1]+gol[jp1+im1]+gol[jm1+i]+gol[jm1+ip1]+gol[jm1+im1];
		newgol[ij] = (1 - ((s>>3 & 0x1) || (s>>2 & 0x1)))  * (s>>1 & 0x1) * ((s&0x1)||gol[ij]);
	}

	for (ij=0; ij<N2; ij++) 
		gol[ij] = newgol[ij];
}

void print (int gol[]) {   /* print the game of life configuration */
	int	ij;
	for (ij=0; ij<N2; ij++) {
		printf ("%c", gol[ij] ? 'x' : ' ');
		if ((ij % N) == N -1) printf ("\n");
	}
}

void initialize (int gol[]) {
	int ij;
	for (ij=0; ij<N2; ij++) {
		gol[ij] = rand() & 0x1;
	}	
}

int main (int argc, char *argv[]) {
	int	gol[N2], i;
    int nsteps = 10000;

    if (argc>1) nsteps = atoi(argv[1]);
	initialize (gol);
    print(gol);
    printf("and after %d steps ......................................................................................\n",nsteps);
	for (i=0; i<nsteps; i++) {
		update (gol);
	}
	print (gol);
    update (gol);
    printf("and after %d steps ......................................................................................\n",nsteps+1);
    print (gol);
}
