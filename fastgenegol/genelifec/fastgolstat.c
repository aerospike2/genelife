//
//  fastgolstat.c
//  fastgenegol
//
//  Created by John McCaskill on 03.09.18.
//  Copyright © 2018 European Center for Living Technology. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>

#define log2N 7

const int N = 0x1 << log2N;
const int N2 = N*N;
const unsigned int Nmask = N - 1;


void update (int gol[], int golstat[]) {
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
		// newgol[ij] = (1 - ((s>>3 & 0x1) | (s>>2 & 0x1)))  * (s>>1 & 0x1) * ((s&0x1)||gol[ij]);
        newgol[ij] = (1 - ((s>>3 & 0x1) | (s>>2 & 0x1)))  & (s>>1 & 0x1) & ((s&0x1)|gol[ij]);
        golstat[(gol[jN+im1]<<7)+(gol[jp1+im1]<<6)+(gol[jp1+i]<<5)+(gol[jp1+ip1]<<4)+(gol[jN+ip1]<<3)+(gol[jm1+ip1]<<2)+(gol[jm1+i]<<1)+gol[jm1+im1]]++;
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

void initialize (int gol[]) {   /* initial GoL field to random values */
	int ij;
	for (ij=0; ij<N2; ij++)  gol[ij] = rand() & 0x1;  /* random starting config */
}

void initializestats (int golstat[]) {   /* initialize the stat counters of local neighbourhoods */
    int i;
    for (i=0; i<256; i++)    golstat[i] = 0;        /* all possible neighborhoods stat counters*/
}

int numones(int x) {
    int i;
    int sum = 0;
    for (i=0; i<8; i++) {
        sum += (x & 0x1);
        x = x >> 1;
    }
    return sum;
}

void printstats (int golstat[]) {
    int i;
    for (i=0; i<256; i++)   {
        printf(" %10x",i);                          /* neighborhoods as hex*/
        if(i%16 == 16-1) printf("\n");
    }
    for (i=0; i<256; i++)   {
        printf(" %10d",numones(i));                 /* number of ones in neighbourhood */
        if(i%16 == 16-1) printf("\n");
    }
    for (i=0; i<256; i++)   {
        if(numones(i)==3) printf(" %10.1f",golstat[i]/900.0);         /* neighborhoods stat counters for 3 live neighbours */
        else printf(" %10s"," ");
        if(i%16 == 16-1) printf("\n");
    }
    for (i=0; i<256; i++)   {
        printf(" %10.1f",golstat[i]/900.0);         /* neighborhoods stat counters*/
        if(i%16 == 16-1) printf("\n");
    }
    
    
}

int main (int argc, char *argv[]) {
	int	gol[N2], i;
    int golstat[256];
    int nsteps = 1000;
    int nafter = 101;    /* reset stat counters at step nafter */

    if (argc>1) nsteps = atoi(argv[1]);
    srand(7654321);
	initialize (gol);            /* random initial pattern */
    printf("initial pattern  .............................................................................................\n");
    print(gol);
	for (i=0; i<nsteps; i++) {   /* nsteps */
        if (i == nafter) initializestats(golstat);
		update (gol, golstat);
	}
    printf("and after %d steps ......................................................................................\n",nsteps);
	print (gol);

    update (gol,golstat);                /* one additional step */
    printf("and after %d steps ......................................................................................\n",nsteps+1);
    print (gol);

    printf("cum stats after %d steps  ...............................................................................\n",nafter);
    printstats (golstat);
    
}
