//
// genelife_size.h
// project genelife
// definition of model size via log2N and derived parameters N,N2,NLM,NLC,Nmask,N2mask: N is length of array side
// declaration of constants specifying array sizes using enum types, to facilitate run time checking of fixed array bounds as well as type checking
//---------------------------------------------------------- copyright -------------------------------------------------------------------------
// Written by John S. McCaskill and Norman H. Packard 2017-2019
//
// First created by John McCaskill on 14.07.2017. Last modified Oct 2019.
// Copyright 2017,2018,2019 European Center for Living Technology. All rights reserved.
//
// This code is distributed in the hope that it will be useful for research purposes, but WITHOUT ANY WARRANTY
// and without even the implied warranty of merchantability or fitness for a particular purpose.
//-----------------------------------------------------------------------------------------------------------------------------------------------
#ifndef genelife_size_h
#define genelife_size_h

#ifndef OEX
#define OEX extern
#endif  /* OEX */

enum {log2N = 7,                	// toroidal array of side length N = 2 to the power of log2N (minimum log2N is 6 i.e. 64x64)
  	  N = 0x1 << log2N,         	// only side lengths powers of 2 allowed to enable efficient implementation of periodic boundaries
 	  N2 = N*N,                 	// number of sites in square-toroidal array
 	  Nmask = N - 1,            	// bit mask for side length, used instead of modulo operation
 	  N2mask = N2 - 1};          	// bit mask for array, used instead of modulo operation
enum {NLM = N2,                     // maximum number of discrete components possible N*N
      NLC = N2<<2};                 // maximum number of connections N*N*4
#endif /* genelife_size_h */
