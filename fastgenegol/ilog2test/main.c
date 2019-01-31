//
//  main.c
//  ilog2test
//
//  Created by John McCaskill on 30/01/2019.
//  Copyright © 2019 ECLT. All rights reserved.
//

#include <stdio.h>
#include <math.h>

extern inline unsigned int mylog2(unsigned int v)// find the log2 of v = power of 2
{
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

extern inline unsigned int mylog2a(unsigned int v)// find the integer log2 of v
{// Adapted by John McCaskill starting from the power of two version at https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
    static const unsigned int b[] = {0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0,
                                     0xFF00FF00, 0xFFFF0000};
    register int k = 4;
    register unsigned int r,s;
    s = ((v & b[k])   != 0); r  = s; v &= s ? b[k] : ~b[k]; r<<=1;
    s = ((v & b[--k]) != 0); r |= s; v &= s ? b[k] : ~b[k]; r<<=1;
    s = ((v & b[--k]) != 0); r |= s; v &= s ? b[k] : ~b[k]; r<<=1;
    s = ((v & b[--k]) != 0); r |= s; v &= s ? b[k] : ~b[k]; r<<=1;
    s = ((v & b[--k]) != 0); r |= s;

    return(r);
}

extern inline unsigned int mylog2b(unsigned int v)// find the integer log2 of v
{// https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
    register unsigned int r; // result of log2(v) will go here
    register unsigned int shift;

    r =     (v > 0xFFFF) << 4; v >>= r;
    shift = (v > 0xFF  ) << 3; v >>= shift; r |= shift;
    shift = (v > 0xF   ) << 2; v >>= shift; r |= shift;
    shift = (v > 0x3   ) << 1; v >>= shift; r |= shift;
                                            r |= (v >> 1);
    return r;
}

extern inline unsigned int mylog2c(unsigned int v)// find the integer log2 of v
{ // https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
    const unsigned int b[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000};
    const unsigned int S[] = {1, 2, 4, 8, 16};
    int i;

    register unsigned int r = 0; // result of log2(v) will go here
    for (i = 4; i >= 0; i--) {   // unroll for speed...
    
        if (v & b[i]) {
            v >>= S[i];
            r |= S[i];
        }
    }
    return r;
}

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
        result = candidateResult*candidateResult <= (n >> shift) ? candidateResult : result;
        shift -= 2;
    }
    return result;
}

int main(int argc, const char * argv[]) {

    unsigned int j;
    unsigned int nside,count=0,count1=0,count2=0,count3=0;
    register unsigned int log2n;

    for(j=0;j<100000;j++)
        for (nside=1;nside<1024;nside++)  {
            for(log2n=0;(1<<log2n)<nside;log2n++);
            count+=log2n;
        }
    
    for(j=0;j<100000;j++)     // 17% of time in profiler  : mine is fastest
        for (nside=1;nside<1024;nside++)  {
            log2n=mylog2a(nside);
            count1+=log2n;
        }
    
    for(j=0;j<100000;j++)     // 25.4 %
        for (nside=1;nside<1024;nside++)  {
            log2n=mylog2b(nside);
            count2+=log2n;
        }

    for(j=0;j<100000;j++)     // 44.5%
        for (nside=1;nside<1024;nside++)  {
            log2n=mylog2c(nside);
            count3+=log2n;
        }

    printf("log sum %d %d %d\n",count1,count2,count3);
    

    for (nside=1;nside<=1024;nside++)  {
        printf("%d %d %d %d\n",nside,mylog2a(nside),mylog2b(nside),mylog2c(nside));
    }

/*
    count=count1=count2=0;
    
    for(j=0;j<100000;j++)
        for (nside=1;nside<1024;nside++)  {
            log2n=(int) sqrt(nside);
            count1+=log2n;
        }
    
    for(j=0;j<100000;j++)
        for (nside=1;nside<1024;nside++)  {
            log2n=integerSqrt(nside);
            count2+=log2n;
        }
    printf("sqrt sum %d %d\n",count1,count2);
    
    for (nside=1;nside<=1000;nside+=7)  {
            printf("%d %d\n",nside,integerSqrt(nside));
    }
    */
    
    return 0;
}
