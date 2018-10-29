2c2
< //  fastgenegol.c
---
> //  fastgenelife.c
25c25
< const int nlog2p0 = 8;              // p0 = 2 to the power of - nlog2p0
---
> const int nlog2p0 = 5;              // p0 = 2 to the power of - nlog2p0
28,30c28,32
< int nsteps = 10000;                 // total number of steps to simulate GoL
< int ndisp = 10000;                  // display GoL every ndisp steps
< int tdisp = 0;                      // extra time delay in ms betwene displays
---
> int nsteps = 100000;                // total number of steps to simulate GoL
> int ndisp  = 10000;                 // display GoL every ndisp steps
> int tdisp  = 0;                     // extra time delay in ms betwene displays
> int rule2mod = 1;                   // whether to modify two live nb rule as well or only three nb rule
> int selection = 1;                  // fitness model: 0 neutral 1 selected gene prob of rule departure 2 presence of replicase gene
44a47
> const long unsigned int m3  = 0x0707070707070707; // for cumcount64c selects counter relevant bits only
51a55,57
> #define CUMCOUNT64C(x, val) {       /* Assumes gene specifies 8 8-bit counters each with max value 7 */  \
>     val = ((x & m3) * h01) >> 56;}  /* left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... */
> 
53,54c59,60
< 	/* update GoL for toroidal field which has side length which is a binary power of 2 */
< 	/* encode without if structures for optimal vector treatment */
---
>     /* update GoL for toroidal field which has side length which is a binary power of 2 */
>     /* encode without if structures for optimal vector treatment */
56c62
< 	int k, nmut, hamming, nb[8], ij, i, j , jp1, jm1, ip1, im1;
---
>     int k, nmut, nones, nones1, nones2, nb[8], ij, i, j , jp1, jm1, ip1, im1;
57a64
>     long unsigned int genef1,genef2,genef3;
63,66c70,73
< 	for (ij=0; ij<N2; ij++) {                                               // loop over all sites of 2D torus with side length N
< 		i = ij & Nmask;  j = ij >> log2N;                                   // row & column
< 		jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
< 		ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                         // toroidal i+1, i-1
---
>     for (ij=0; ij<N2; ij++) {                                               // loop over all sites of 2D torus with side length N
> 	i = ij & Nmask;  j = ij >> log2N;                                   // row & column
> 	jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
> 	ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                         // toroidal i+1, i-1
69c76
< 		s = gol[nb[0]]+gol[nb[1]]+gol[nb[2]]+gol[nb[3]]+gol[nb[4]]+gol[nb[5]]+gol[nb[6]]+gol[nb[7]]; // number of live nbs
---
> 	s = gol[nb[0]]+gol[nb[1]]+gol[nb[2]]+gol[nb[3]]+gol[nb[4]]+gol[nb[5]]+gol[nb[6]]+gol[nb[7]]; // number of live nbs
71,95c78,122
<         genediff = s2or3 * (golg[nb[nb1i&0x7]]^golg[nb[(nb1i>>4)&0x7]]^golg[nb[(nb1i>>8)&0x7]]);  // gene difference sequence based on xor for 2 or 3 live nbs
<         POPCOUNT64C(genediff, hamming);
<         
<         // compute random events r1 and r2, neighbor selection ng and single mutation position nmut from same 64-bit random number
<         RAND128P(randnr);                                                   // expansion inline so compiler recognizes auto-vectorization options
<         nlog2p = nlog2p0 + hamming;                                         // real factor alpha not possible here, could do integer conversion
<         pmask = (0x1<<nlog2p) - 1;                                          // probability mask for deviation from gol rules given local hamming
<         randnr1 = randnr & pmask;                                           // extract bits from randnr for random trial for 0 on pmask
<         r1 = ((randnr1 | ((randnr1^pmask)+1)) >> nlog2p)&0x1;               // 1 if lowest nlog2p bits of randnr zero, else zero : i.e. 1 with chance 1/2^nlog2p
<         randnr2 = (randnr >> 16) & pmutmask;                                // extract bits from randnr for random trial for 0 on pmask
<         r2 = ((randnr2 | ((randnr2^pmutmask)+1)) >> nlog2pmut)&0x1;         // 1 if lowest nlog2pmut bits of randnr zero, else zero
<         nmut = (randnr >> 32) & 0x3f;                                       // choose mutation position for length 64 gene
<         ng = (randnr >> 38) & 0x3;                                          // 0, 1, 2 or 3 with probs each 1/4 : next 4 lines converts this 0,1,2 with prob 1/3
<         r3 = ((ng+1) >> 2) & 0x1;                                           // 1 if ng == 3 (invalid value) otherwise zero : prob. is 1/4
<         ngx += r3;                                                          // increment external ng counter on such exceptions
<         ngx = (1-(((ngx+1)>> 2)&0x1))*ngx;                                  // modulo 3 counter without division for ng==3 exceptions
<         ng = (s&1)*((1-r3)*ng+r3*ngx)+(1-(s&1))*(ng&1);                     // use counter value mod 3 if exception, else random 0,1,2 for s==3; random 0,1 for 2
<         
<                                   // number of 1s in genediff is Hamming distance
< 
<         newgene = golg[nb[(nb1i>>(ng<<2))& 7]];                             // pick new gene as one of three neighbor live site genes
<         newgene = newgene ^ (r2*(0x1L<<nmut));                              // introduce single mutation with probability pmut = probmut
<         birth = (1-gol[ij])&((s&1)^r1)&0x1;                                 // assuming 2or3 live nbs, birth (value 1) if empty and (s==3 xor r1)
<         newgol[ij]  = s2or3 * ( gol[ij] | birth );                          // new game of life cell value
<         newgolg[ij] = s2or3 * ( gol[ij]*golg[ij]+birth*newgene);            // dies if not 2or3, else old if alive, else new gene if 3 nbs
---
>         if (s2or3 == 1) {                                                   // if 2 or 3 neighbors alive
>             RAND128P(randnr);                                                 // expansion inline so compiler recognizes auto-vectorization options
>                                                                               // compute random neighbor selection ng from one 64-bit random number
>             ng = (randnr >> 54) & 0x3;                                          // 0, 1, 2 or 3 with probs each 1/4 : next 4 lines converts this 0,1,2 with prob 1/3
>             r3 = (ng == 3 ? 1: 0);                                              // 1 if ng == 3 (invalid value) otherwise zero : prob. is 1/4
>             ngx += r3;                                                          // increment external ng counter on such exceptions
>             ngx = (ngx == 3 ? 0 : ngx);                                         // modulo 3 counter without division for ng==3 exceptions
>             ng = (s&1)?(r3?ngx:ng):(ng&1);                                      // for s==3 use counter value mod 3 if exception, else rand 0,1,2; for s=2 rand 0,1
>             newgene = golg[nb[(nb1i>>(ng<<2))& 7]];                             // pick new gene as one of three neighbor
> 	    // compute nones for different selection models
>             if (selection == 0) {                                               // neutral model : GoL rule departures depend only on seq diversity
>                 genediff = golg[nb[nb1i&0x7]]^golg[nb[(nb1i>>4)&0x7]]^golg[nb[(nb1i>>8)&0x7]];  // gene difference seq based on xor
>                 POPCOUNT64C(genediff, nones);}                                    // number of 1s in genediff is Hamming distance
>             else if (selection == 1) {                                          // non-neutral model with selection for rule departure probability
>                 genediff = newgene;                                               // place the new gene in genediff for count of number of ones
>                 POPCOUNT64C(genediff, nones);                                     // number of ones in new gene determines fitness
>                 nones = (nones < 16) ? 0 : (nones - 23);}                           // 0 if < 16 otherwise nones-23
>             else {                                                              // non-neutral model based on presence of replicase gene
>                 genef1 = golg[nb[nb1i&0x7]];                                      // gene difference seq based on xor
>                 CUMCOUNT64C(genef1, nones);                                       // number of ones in new gene determines fitness
>                 genef2 = golg[nb[(nb1i>>4)&0x7]];
>                 CUMCOUNT64C(genef2, nones1);                                      // number of ones in new gene determines fitness
>                 genef3 = golg[nb[(nb1i>>8)&0x7]];
>                 CUMCOUNT64C(genef3, nones2);                                      // number of ones in new gene determines fitness
>                 nones = nones < nones1 ? nones : nones1;
>                 nones = nones < nones2 ? nones : nones2;
>                 nones = nones < 16 ? 0 : (nones - 23);}                           // 0 if < 16 otherwise nones-23
> 	    // compute departure and mutation events, mutation position nmut
>             nlog2p = nlog2p0 + nones;                                           // real factor alpha not possible here, could do integer conversion
>             pmask = (0x1<<nlog2p) - 1L;                                          // probability mask for deviation from gol rules given local hamming
>             randnr1 = randnr & pmask;                                           // extract bits from randnr for random trial for 0 on pmask
>             r1 = randnr1?0L:1L;                                                   // 1 if lowest nlog2p bits of randnr zero, else zero : i.e. 1 with chance 1/2^nlog2p
>             randnr2 = (randnr >> 24) & pmutmask;                                // extract bits from randnr for random trial for 0 on pmask
>             r2 = randnr2?0:1;                                                   // 1 if lowest nlog2pmut bits of randnr zero, else zero
>             nmut = (randnr >> 48) & 0x3f;                                       // choose mutation position for length 64 gene : from bits 32:37 of randnr
> 	    // complete calculation of newgol and newgolg, including mutation
>             newgene = newgene ^ (r2*(0x1L<<nmut));                              // introduce single mutation with probability pmut = probmut
>             birth = (0x1L-gol[ij])&((s&1L)^(r1&rule2mod))&0x1;                                 // assuming 2or3 live nbs, birth (value 1) if empty and (s==3 xor r1)
>             newgol[ij]  =  gol[ij] | birth ;                                    // new game of life cell value
>             newgolg[ij] =  gol[ij]*golg[ij]+birth*newgene;}                     // dies if not 2or3, else old if alive, else new gene if 3 nbs
>         else {                                                              // else not 2 or 3 live neighbors, 0 values
>             newgol[ij]  = 0;                                                    // new game of life cell value
>             newgolg[ij] = 0;}                                                   // dies if not 2or3
>         emptysites = emptysites + newgol[ij];                               // keep track of empty sites, same information as total activity of occupied sites
>     }
97c124,126
<         emptysites = emptysites + newgol[ij];
---
>     for (ij=0; ij<N2; ij++) {
> 	gol[ij] = newgol[ij];        // copy new gol config to old one
> 	golg[ij] = newgolg[ij];      // copy new genes to old genes
98a128
> }
100,102c130,193
< 	for (ij=0; ij<N2; ij++) {
<         gol[ij] = newgol[ij];        // copy new gol config to old one
<         golg[ij] = newgolg[ij];      // copy new genes to old genes
---
> // like update, but has plane, newplane (for both gol, golg) both as args.
> void update2 (long unsigned int gol[], long unsigned int newgol[], long unsigned int golg[], long unsigned int newgolg[]) {
>     /* update GoL for toroidal field which has side length which is a binary power of 2 */
>     /* encode without if structures for optimal vector treatment */
> 
>     int k, nmut, nones, nones1, nones2, nb[8], ij, i, j , jp1, jm1, ip1, im1;
>     long unsigned int s, s2or3, nb1i, randnr, randnr1, randnr2, ng, r1, r2, r3, nlog2p, pmask, genediff, birth, newgene;
>     long unsigned int genef1,genef2,genef3;
>     static long unsigned int  pmutmask = (0x1 << nlog2pmut) - 1;
>     static long unsigned int ngx = 0;
> 
> 
>     for (ij=0; ij<N2; ij++) {                                               // loop over all sites of 2D torus with side length N
> 	i = ij & Nmask;  j = ij >> log2N;                                   // row & column
> 	jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
> 	ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                         // toroidal i+1, i-1
>         nb[0]=j*N+ip1; nb[1]=j*N+im1; nb[2]=jp1+i; nb[3]=jp1+ip1; nb[4]=jp1+im1; nb[5]=jm1+i; nb[6]=jm1+ip1; nb[7]=jm1+im1; //nbs
>         for (k=0,nb1i=0;k<8;k++) nb1i = (nb1i << (gol[nb[k]]<<2)) + (gol[nb[k]]*k);   // packs non-zero nb indices in first up to 8*4 bits
> 	s = gol[nb[0]]+gol[nb[1]]+gol[nb[2]]+gol[nb[3]]+gol[nb[4]]+gol[nb[5]]+gol[nb[6]]+gol[nb[7]]; // number of live nbs
>         s2or3 = (1 - (((s>>3)&1) | ((s>>2)&1))) * (s>>1 & 1);               // 1 if 2 or 3 neighbors are alive
>         if (s2or3 == 1) {                                                   // if 2 or 3 neighbors alive
>             RAND128P(randnr);                                                 // expansion inline so compiler recognizes auto-vectorization options
>                                                                               // compute random neighbor selection ng from one 64-bit random number
>             ng = (randnr >> 54) & 0x3;                                          // 0, 1, 2 or 3 with probs each 1/4 : next 4 lines converts this 0,1,2 with prob 1/3
>             r3 = (ng == 3 ? 1: 0);                                              // 1 if ng == 3 (invalid value) otherwise zero : prob. is 1/4
>             ngx += r3;                                                          // increment external ng counter on such exceptions
>             ngx = (ngx == 3 ? 0 : ngx);                                         // modulo 3 counter without division for ng==3 exceptions
>             ng = (s&1)?(r3?ngx:ng):(ng&1);                                      // for s==3 use counter value mod 3 if exception, else rand 0,1,2; for s=2 rand 0,1
>             newgene = golg[nb[(nb1i>>(ng<<2))& 7]];                             // pick new gene as one of three neighbor
> 	    // compute nones for different selection models
>             if (selection == 0) {                                               // neutral model : GoL rule departures depend only on seq diversity
>                 genediff = golg[nb[nb1i&0x7]]^golg[nb[(nb1i>>4)&0x7]]^golg[nb[(nb1i>>8)&0x7]];  // gene difference seq based on xor
>                 POPCOUNT64C(genediff, nones);}                                    // number of 1s in genediff is Hamming distance
>             else if (selection == 1) {                                          // non-neutral model with selection for rule departure probability
>                 genediff = newgene;                                               // place the new gene in genediff for count of number of ones
>                 POPCOUNT64C(genediff, nones);                                     // number of ones in new gene determines fitness
>                 nones = (nones < 16) ? 0 : (nones - 23);}                           // 0 if < 16 otherwise nones-23
>             else {                                                              // non-neutral model based on presence of replicase gene
>                 genef1 = golg[nb[nb1i&0x7]];                                      // gene difference seq based on xor
>                 CUMCOUNT64C(genef1, nones);                                       // number of ones in new gene determines fitness
>                 genef2 = golg[nb[(nb1i>>4)&0x7]];
>                 CUMCOUNT64C(genef2, nones1);                                      // number of ones in new gene determines fitness
>                 genef3 = golg[nb[(nb1i>>8)&0x7]];
>                 CUMCOUNT64C(genef3, nones2);                                      // number of ones in new gene determines fitness
>                 nones = nones < nones1 ? nones : nones1;
>                 nones = nones < nones2 ? nones : nones2;
>                 nones = nones < 16 ? 0 : (nones - 23);}                           // 0 if < 16 otherwise nones-23
> 	    // compute departure and mutation events, mutation position nmut
>             nlog2p = nlog2p0 + nones;                                           // real factor alpha not possible here, could do integer conversion
>             pmask = (0x1<<nlog2p) - 1L;                                          // probability mask for deviation from gol rules given local hamming
>             randnr1 = randnr & pmask;                                           // extract bits from randnr for random trial for 0 on pmask
>             r1 = randnr1?0L:1L;                                                   // 1 if lowest nlog2p bits of randnr zero, else zero : i.e. 1 with chance 1/2^nlog2p
>             randnr2 = (randnr >> 24) & pmutmask;                                // extract bits from randnr for random trial for 0 on pmask
>             r2 = randnr2?0:1;                                                   // 1 if lowest nlog2pmut bits of randnr zero, else zero
>             nmut = (randnr >> 48) & 0x3f;                                       // choose mutation position for length 64 gene : from bits 32:37 of randnr
> 	    // complete calculation of newgol and newgolg, including mutation
>             newgene = newgene ^ (r2*(0x1L<<nmut));                              // introduce single mutation with probability pmut = probmut
>             birth = (0x1L-gol[ij])&((s&1L)^(r1&rule2mod))&0x1;                                 // assuming 2or3 live nbs, birth (value 1) if empty and (s==3 xor r1)
>             newgol[ij]  =  gol[ij] | birth ;                                    // new game of life cell value
>             newgolg[ij] =  gol[ij]*golg[ij]+birth*newgene;}                     // dies if not 2or3, else old if alive, else new gene if 3 nbs
>         else {                                                              // else not 2 or 3 live neighbors, 0 values
>             newgol[ij]  = 0;                                                    // new game of life cell value
>             newgolg[ij] = 0;}                                                   // dies if not 2or3
>         emptysites = emptysites + newgol[ij];                               // keep track of empty sites, same information as total activity of occupied sites
165c256
<     int ij, k, ijlast, nspecies, counts[N2];
---
>     int ij, k, ijlast, nspecies, counts[N2], nones, fitness;
191c282,293
<         printf("count species %d with gene %lx has counts %lu\n",k, golgsc[k][0],golgsc[k][1]);
---
>         last = golgsc[k][0];
>         
>         if (selection == 0) {                                               // neutral model : GoL rule departures depend only on seq diversity
>                 POPCOUNT64C(last, nones);
>                 fitness = nlog2p0;}
>         else if (selection == 1) {                                          // non-neutral model with selection for rule departure probability
>                 POPCOUNT64C(last, nones);                                     // number of ones in new gene determines fitness
>                 fitness = nlog2p0 + ((nones < 16) ? 0 : (nones - 23));}       // 0 if < 16 otherwise nones-23
>         else {                                                              // non-neutral model based on presence of replicase gene
>                 POPCOUNT64C(last, nones);                                     // number of ones in new gene determines fitness
>                 fitness = nlog2p0 + ((nones < 16) ? 0 : (nones - 23));}       // 0 if < 16 otherwise nones-23
>         printf("count species %d with gene %lx has counts %lu and %d ones, fitness %d\n",k, golgsc[k][0],golgsc[k][1],nones,fitness);
