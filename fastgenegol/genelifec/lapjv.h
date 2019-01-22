#ifndef LAPJV_H
#define LAPJV_H

#define LARGE 1000000

#if !defined TRUE
#define TRUE 1
#endif
#if !defined FALSE
#define FALSE 0
#endif

#define NEW(x, t, n) if ((x = (t *)malloc(sizeof(t) * (n))) == 0) { return -1; }
#define FREE(x) if (x != 0) { free(x); x = 0; }
#define SWAP_INDICES(a, b) { int_t _temp_index = a; a = b; b = _temp_index; }

#if 0
#include <assert.h>
// #define LAPMOD_ASSERT(cond) assert(cond)
#define LAPMOD_ASSERT(cond,nr,a,b) { if (!(cond)) fprintf(stderr,"LAPMOD_ASSERT %d fails for %d,%d \n",(nr),(a),(b)); assert(cond); }     /* modified to check asserts */
#define PRINTF(fmt, ...) printf(fmt, ##__VA_ARGS__)
#define PRINT_COST_ARRAY(a, n) \
    while (1) { \
        fprintf(stderr,#a" = ["); \
        if ((n) > 0) { \
            fprintf(stderr,"%f", (a)[0]); \
            for (uint_t j = 1; j < n; j++) { \
                fprintf(stderr,", %f", (a)[j]); \
            } \
        } \
        fprintf(stderr,"]\n"); \
        break; \
    }
#define PRINT_INDEX_ARRAY(a, n) \
    while (1) { \
        fprintf(stderr,#a" = ["); \
        if ((n) > 0) { \
            fprintf(stderr,"%d", (a)[0]); \
            for (uint_t j = 1; j < n; j++) { \
                fprintf(stderr,", %d", (a)[j]); \
            } \
        } \
        fprintf(stderr,"]\n"); \
        break; \
    }
#else
#include <assert.h>                   // added to check asserts
// #define LAPMOD_ASSERT(cond)
#define LAPMOD_ASSERT(cond,nr,a,b) { if (!(cond)) fprintf(stderr,"LAPMOD_ASSERT %d fails for %d,%d \n",(nr),(a),(b)); assert(cond); }     /* modified to check asserts */
#define PRINTF(fmt, ...)
#define PRINT_COST_ARRAY(a, n)
#define PRINT_INDEX_ARRAY(a, n)
#endif


typedef int int_t;
typedef unsigned int uint_t;
typedef double cost_t;
typedef char boolean;
typedef enum fp_t { FP_1 = 1, FP_2 = 2, FP_DYNAMIC = 3 } fp_t;

/* extern int_t lapjv_internal(
    const uint_t n, cost_t *cost[],
    int_t *x, int_t *y); */          // not included or used in genelife, always sparse

extern int_t lapmod_internal(
    const uint_t n, cost_t *cc, uint_t *ii, uint_t *kk,
    int_t *x, int_t *y, fp_t fp_version);

#endif // LAPJV_H
/*
Copyright (c) 2012-2017, Tomas Kazmar
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
