//
//  compare.c
//  project genelife
//
//  custom comparison functions, primarily for use in sorting
//---------------------------------------------------------- copyright --------------------------------------------------------------------------------
//  Written by John S. McCaskill and Norman H. Packard 2017-2019
//
//  First created by John McCaskill on 14.07.2017. Last modified Oct 2019.
//  Copyright 2017,2018,2019 European Center for Living Technology. All rights reserved.
//
//  This code is distributed in the hope that it will be useful for research purposes, but WITHOUT ANY WARRANTY
//  and without even the implied warranty of merchantability or fitness for a particular purpose.
//------------------------------------------------------------------------------------------------------------------------------------------------------
//
#include "genelife.h"
//
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
