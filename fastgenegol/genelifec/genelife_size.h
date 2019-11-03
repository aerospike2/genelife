enum {log2N = 7,                	// toroidal array of side length N = 2 to the power of log2N (minimum log2N is 6 i.e. 64x64)
  	  N = 0x1 << log2N,         	// only side lengths powers of 2 allowed to enable efficient implementation of periodic boundaries
 	  N2 = N*N,                 	// number of sites in square-toroidal array
 	  Nmask = N - 1,            	// bit mask for side length, used instead of modulo operation
 	  N2mask = N2 - 1};          	// bit mask for array, used instead of modulo operation
enum {NLM};

