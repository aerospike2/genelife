//
//  lsb64.c
//  fastgenegol
//
//  Created by John McCaskill on 27/11/2018.
//  Copyright © 2018 Ruhr Universität Bochum. All rights reserved.
//

#include <stdio.h>
#include <stdint.h>

int main () {
unsigned int i,j,cnt,cnt2;
uint64_t u,v,val,sum;

#define FIRST1INDEX(v, c) {                    /* starting point 64bit from Sean Eron Anderson https://graphics.stanford.edu/~seander/bithacks.html#ZerosOnRightParallel */  \
    uint64_t mmmm,mmmq;                        /* note also arguments must be of types uint64_t and int respectivley */ \
    int cccc;                                  /* takes on successive integer values 32,16,84,2,1 */ \
    int t1f0;                                  /* logical to integer variable true=one false=zero : if a 1 in v under mask mmmm */ \
    c=64;                                      /* contains count of number of zeros on right of last one */ \
    mmmm=~0ull;                                /* all ones, mask from previous stage */ \
    for (cccc=c>>1;cccc>0;cccc>>=1) {          /* loop over cccc goes 32,16,8,4,2,1 */ \
        mmmq = mmmm;                           /* mmmq is to be the mask used to query if a one is under it at this stage, start with old mask */ \
        mmmq &= mmmm^(mmmm<<cccc);             /* divided high part of mask into two equal parts, taking lower part */ \
        t1f0 = v&mmmq?1:0;                     /* one if a one under the query mask, zero otherwise */ \
        mmmm=mmmq^((1-t1f0)*mmmm);             /* the new mask for next stage is the query mask if a one is under it, otherwise the other half of mmmm */ \
        c-=t1f0*cccc;                          /* the right zero counter is decremented by the length of the current interval cccc if a one was under mask */ \
                                               /* fprintf(stderr,"in FIRST1INDEX at cccc = %d with mmmm = %llx and c = %d\n",cccc,mmmm,c); */ \
    }                                          /* note that Anderson's algorithm was incorrect */ \
    if (v) c--;                                /* only the case with no ones at all avoids the additional decrement caused by the presence of a one */ \
}

#define FIRST1INDEX2(v, c) {                   /* starting point 64bit from Sean Eron Anderson https://graphics.stanford.edu/~seander/bithacks.html#ZerosOnRightParallel */  \
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
    }                                          /* note that Anderson's algorithm was incorrect */ \
}
u=42ull;
for(j=sum=0;j<10000;j++) {
  u++;
  for (v=u,i=0;i<100000;i++) {
    val = v*i*11400714819323198549ull;
    // FIRST1INDEX(val, cnt)
    FIRST1INDEX2(val, cnt)
    // if(cnt!=cnt2) printf("error at step j %d i %d\n",j,i);
    sum+=cnt;
  }
}

printf("sum of values %llx\n",sum);
}
