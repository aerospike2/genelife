//
//  quadhash.c
//  fastgenegol
//
//  Created by John McCaskill on 17/12/2018.
//  Copyright © 2018 Ruhr Universität Bochum. All rights reserved.
//

#include "quadhash.h"
const int N = 512;
const int N2 = N*N;
const int Nmask =N-1;
const int N2mask =N2-1;
//.......................................................................................................................................................
#define deltaxyl(ij,l,x,y)  (ij - (ij&(Nmask>>l)) + (((ij+(x))&(Nmask>>l)) + (y)*(N>>l))) & (N2mask>>(l<<1))

extern inline void quadimage(uint64_t gol[],uint64_t golp[]) {              // routine to pack 8x8 subarrays into single words
    unsigned int ij,ij8,k;
    
    for(k=0;k<64;k++)
        for (ij=ij8=0;ij<N2;ij+=8,ij8++) {
             golp[ij8] |= gol[deltaxyl(ij,k&0x7,k>>3,0)]<<k;                         // 8x8 packed arrays
             ij+= (ij+1)&((N<<3)-1) ? 0 : (8-1)*N;
        }
}

struct node {
   struct node *next;                                                             // hash-linked list of elements with same key
   struct node *nw, *ne, *sw, *se;                                                // constant pointers; nw != 0 for non terminal nodes
   struct node *res;                                                              // store cached value associated with the node
};
typedef struct node node;

struct leaf {
   node *next;                                                             // ptr to next elt in hash-linked list of elts with same key
   node *isnode;                                                           // isnode = 0 for leaves
   uint64_t pattern;                                                       // store 64 bit pattern associated with the leaf
   uint64_t reserve:                                                       // reserve value, pads to 4 words, useful for tracing etc
};
typedef struct leaf leaf;

uint64_t hashmask;

node **hashtab;     // hash table of nodes : pointers to lists of chained nodes that all map to the same hash addresses
node **stack;       // stack of all the roots
node *freenodes;    // list of free nodes (allocated in sets of 1000
node * nodeblocks;  //
int gsp,maxmem;
int stacksize;
int okaytogc ;
uint64_t alloced,totalthings;

typedef uint64_t * uint64p_t;

extern inline uint64_t node_hash(const void *a, const void *b, const void *c, const void *d) {
   uint64_t r = (65537*(uint64_t)(d)+257*(uint64_t)(c)+17*(uint64_t)(b)+5*(uint64_t)(a));
   r += (r >> 11);
   return r ;
}

node * save(node *n) {  // This routine marks a node as needing to be saved.
   if (gsp >= stacksize) {  // stacks needs enlarging
      int nstacksize = stacksize * 2 + 100;
      alloced += sizeof(node *)*(nstacksize-stacksize);
      stack = (node **)realloc(stack, nstacksize * sizeof(node *));
      if (stack == 0) {
        fprintf(stderr,"Run out of memory\n");
        exit(1);
      }
      stacksize = nstacksize;
   }
   stack[gsp++] = n ;
   return n ;
}

#define HASHMOD(a) ((a)&(hashmask))

static uint64_t nexthashsize(uint64_t i) {
   while ((i & (i - 1))) // stops when number is next power of two ie in the form 1000...0
      i += (i & (1 + ~i)) ;
   return i ;
}

node *newnode() {                // free nodes kept in a linked list, allocated 1000 at a time
   node *r ;
   if (freenodes == 0) {
      int i ;
      freenodes = (node *) calloc(1001, sizeof(node)) ;
      if (freenodes == 0) {
         fprintf(stderr,"Run out of memory; try reducing the hash memory limit.\n") ;
         exit(1);
      }
      alloced += 1001 * sizeof(node) ;
      freenodes->next = nodeblocks ;
      nodeblocks = freenodes++ ;
      for (i=0; i<999; i++) {
         freenodes[1].next = freenodes ;
         freenodes++ ;
      }
      totalthings += 1000 ;
   }
   if (freenodes->next == 0 && alloced + 1000 * sizeof(node) > maxmem && okaytogc) do_gc(0) ;

   r = freenodes ;
   freenodes = freenodes->next ;
   return r ;
}

node *find_node(node *nw, node *ne, node *sw, node *se) {
   node *p ;
   uint64_t h;
   
   h = node_hash(nw,ne,sw,se);
   node *pred = 0;
   h &= hashmask;
   for (p=hashtab[h]; p; p = p->next) {
      if (nw == p->nw && ne == p->ne && sw == p->sw && se == p->se) {
         if (pred) { /* move this one to the front */
            pred->next = p->next;
            p->next = hashtab[h];
            hashtab[h] = p;
         }
         return save(p);
      }
      pred = p;
   }
   p = newnode();
   p->nw = nw;
   p->ne = ne;
   p->sw = sw;
   p->se = se;
   p->value = 0;
   p->next = hashtab[h];
   hashtab[h] = p;
   hashpop++;
   save(p);
   if (hashpop > hashlimit)
      resize();
   return p;
}

leaf *find_leaf(unsigned short nw, unsigned short ne,
                                  unsigned short sw, unsigned short se) {
   leaf *p ;
   leaf *pred = 0 ;
   g_uintptr_t h = leaf_hash(nw, ne, sw, se) ;
   h = HASHMOD(h) ;
   for (p=(leaf *)hashtab[h]; p; p = (leaf *)p->next) {
      if (nw == p->nw && ne == p->ne && sw == p->sw && se == p->se &&
          !is_node(p)) {
         if (pred) {
            pred->next = p->next ;
            p->next = hashtab[h] ;
            hashtab[h] = (node *)p ;
         }
         return (leaf *)save((node *)p) ;
      }
      pred = p ;
   }
   p = newleaf() ;
   p->nw = nw ;
   p->ne = ne ;
   p->sw = sw ;
   p->se = se ;
   leafres(p) ;
   p->isnode = 0 ;
   p->next = hashtab[h] ;
   hashtab[h] = (node *)p ;
   hashpop++ ;
   save((node *)p) ;
   if (hashpop > hashlimit)
      resize() ;
   return p ;
}

void hlifealgo() {
   int i ;
/*
 *   The population of one-bits in an integer is one more than the
 *   population of one-bits in the integer with one fewer bit set,
 *   and we can turn off a bit by anding an integer with the next
 *   lower integer.
 */
   if (shortpop[1] == 0)
      for (i=1; i<65536; i++)
         shortpop[i] = shortpop[i & (i - 1)] + 1 ;
   hashprime = nexthashsize(1000) ;
#ifndef PRIMEMOD
   hashmask = hashprime - 1 ;
#endif
   hashlimit = (g_uintptr_t)(maxloadfactor * hashprime) ;
   hashpop = 0 ;
   hashtab = (node **)calloc(hashprime, sizeof(node *)) ;
   if (hashtab == 0)
     lifefatal("Out of memory (1).") ;
   alloced = hashprime * sizeof(node *) ;
   ngens = 0 ;
   stacksize = 0 ;
   halvesdone = 0 ;
   nzeros = 0 ;
   stack = 0 ;
   gsp = 0 ;
   maxmem = 256 * 1024 * 1024 ;
   freenodes = 0 ;
   okaytogc = 0 ;
   totalthings = 0 ;
   nodeblocks = 0 ;
   zeronodea = 0 ;
   ruletable = hliferules.rule0 ;
/*
 *   We initialize our universe to be a 16-square.  We are in drawing
 *   mode at this point.
 */
   root = (node *)newclearednode() ;
   population = 0 ;
   generation = 0 ;
   increment = 1 ;
   setincrement = 1 ;
   nonpow2 = 1 ;
   pow2step = 1 ;
   llsize = 0 ;
   depth = 3 ;
   hashed = 0 ;
   popValid = 0 ;
   needPop = 0 ;
   inGC = 0 ;
   cacheinvalid = 0 ;
   gccount = 0 ;
   gcstep = 0 ;
   running_hperf.clear() ;
   inc_hperf = running_hperf ;
   step_hperf = running_hperf ;
   softinterrupt = 0 ;
}



//.......................................................................................................................................................
extern inline int comparequad(struct node * golq1, struct node * golq2) {         // recursive comparison of quadtrees at same level of hierarchy

    if (golq1->nw && golq2->nw)
        return(comparequad(golq1->nw,golq2->nw) && comparequad(golq1->ne,golq2->ne) &&
               comparequad(golq1->sw,golq2->sw) && comparequad(golq1->se,golq2->se));
    else if (golq1->nw || golq2->nw) return(0);
    else if (golq1->value == golq2->value)
        return(1);
    else return(0);
}

//.......................................................................................................................................................


quadnode * quadimage(uint64_t gol[]) {                                          // routine to generate a quadtree for an entire binary image of long words
                                                                                // makes use of global linear and quadratic size variables N and N2 for gol
    unsigned int ij,ij1,n;
    uint64_t golp[N2>>6];
    quadnode * golq[N2>>8];
    
    pack64neighbors(gol,golp);

    n=N>>3;
    for (ij=ij1=0;ij<n*n;ij+=2,ij1++) {
        golq[ij1]=hash_patt_find(golp[(ij+n)],golp[(ij+n)+1],golp[ij],golp[ij+1]); // hash_patt_find(nw,ne,sw,se) adds quad leaf entry if pattern not found
        ij+= ((ij+2)&(n-1)) ? 0 : n;                                            // skip odd rows since these are the northern parts of quads generated on even rows
    }
    
    for(n >>= 1; n>1; n>>= 1) {
        for (ij=ij1=0;ij<n*n;ij+=2,ij1++) {
            golq[ij1]=hash_node_find(golq[(ij+n)],golq[(ij+n)+1],golq[ij],golq[ij+1]); // hash_node_find(nw,ne,sw,se) adds quad node entry if node not found
            ij+= ((ij+2)&(n-1)) ? 0 : n;                                        // skip odd rows since these are the northern parts of quads generated on even rows
        }
    }
    if(golq[0]!=NULL)
        if(golq[0]->hits > 1) fprintf(stderr,"image already found at t = %d\n",golq[0]->firsttime);
    return(golq[0]);
}
