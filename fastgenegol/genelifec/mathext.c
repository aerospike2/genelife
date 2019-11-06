//
//  mathext.c
//  project genelife
//
//  some simple auxiliary math functions which may not be defined by the system
//---------------------------------------------------------- copyright --------------------------------------------------------------------------------
//  Written by John S. McCaskill and Norman H. Packard 2017-2019
//
//  First created by John McCaskill on 14.07.2017. Last modified Oct 2019.
//  Copyright 2017,2018,2019 European Center for Living Technology. All rights reserved.
//
//  This code is distributed in the hope that it will be useful for research purposes, but WITHOUT ANY WARRANTY
//  and without even the implied warranty of merchantability or fitness for a particular purpose.
//------------------------------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#define INLINE
//
// integerSqrt          direct bit processing algorithm to implement integer sqrt (largest integer smaller than sqrt) : but floating point sqrt is faster
// log2r                fast integer logarithm working only for arguments which are powers of 2 (not used but slightly faster than log2a when applicable)
// log2lower            fast integer logarithm working for all integers : largest integer smaller than or equal to logarithm base 2 of argument
// log2upper            fast integer logarithm working for all integers : smallest integer larger than or equal to logarithm base 2 of argument
// sqrtupper            fast integer sqrt working for all integers : smallest integer larger than or equal to sqrt of argument
// randprob             random event with probability determined by a 32 bit unsigned integer iprob as iprob / 2^32 using RAND128 and uint64_t
//--------------------------------------------------------------- mathematical fns ----------------------------------------------------------------------
extern INLINE int integerSqrt(int n) {                  // the largest integer smaller than the square root of n (n>=0)
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
extern INLINE unsigned int log2r(unsigned int v) { // find the log2 of v = power of 2, warning wrong answers for v not power of 2
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
extern INLINE unsigned int log2lower(unsigned int v) { // find the integer log2 of v (lower): works for all 32 bit integer v > 0
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
extern INLINE unsigned int log2upper(unsigned int v) { // find the integer log2 of v (upper): works for all 32 bit integer v > 0

    return v ? 1 + log2lower(v-1) : 0;              // assuming as in this implementation that log2lower(0) == 0
}
//.......................................................................................................................................................
extern INLINE unsigned int sqrtupper(unsigned int v) { // smallest integer larger than or equal to sqrt of argument

    return v ? 1 + (int) sqrt(v-1) : 0;
}
//..............................................................  randprob  .............................................................................
extern INLINE uint64_t randprob(unsigned int uprob, unsigned int randnr) {
    return(randnr < uprob ? 1ull : 0ull);
}
