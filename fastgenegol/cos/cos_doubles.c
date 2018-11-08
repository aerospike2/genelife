//
//  cos_doubles.c
//  fastgenegol
//
//  Created by John McCaskill on 24.07.17.
//  Copyright © 2017 Ruhr Universiät Bochum. All rights reserved.
//

#include "cos_doubles.h"

#include <math.h>

/*  Compute the cosine of each element in in_array, storing the result in
 *  out_array. */
void cos_doubles(double * in_array, double * out_array, int size){
    int i;
    for(i=0;i<size;i++){
        out_array[i] = cos(in_array[i]);
    }
}
