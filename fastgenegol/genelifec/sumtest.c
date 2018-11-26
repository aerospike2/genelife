//
//  sumtest.c
//  fastgenegol
//
//  Created by John McCaskill on 26/11/2018.
//  Copyright © 2018 Ruhr Universität Bochum. All rights reserved.
//

//  code snippet to test update_gol64 3D sum formation

    unsigned int sum0,sum1,sum2,sum3,j,j1;  // debugging only
    uint64_t sumsdebug0[4],sumsdebug1[4],sumsdebug2[4],sumsdebug3[4];
         //***********************************************************************
        // check s6
        if(rulemod && b3d) {
            for(j=0;j<64;j++) {
                k=j&0x3;
                j1=j>>2;
                sum0 = (s6>>j)&1ull;
                sum1 = (sums3D[k]>>(j1<<2))&0xf;
                if(sum0 && (sum1!=6)) fprintf(stderr,"check failed: s6 case for which sum is not 6 but %d\n",sum1);
                plane=j;
                sum1=(sums16[plane&0x3]>>(plane&0x3c))&0xf;
                plane = (j+1) & 0x3f;
                sum1+=sum2=(sums16[plane&0x3]>>(plane&0x3c))&0xf;
                plane = (j-1) & 0x3f;
                sum1+=sum3=(sums16[plane&0x3]>>(plane&0x3c))&0xf;
                if ((b3d>>j)&1ull)   // central state should be zero (yes)
                    if(sum0 && (sum1!=6)) fprintf(stderr,"check failed: s6 case j %d for which plane sum is not 6 from 3D %d %d %d (%d %d,%d,%d) but %d (%d,%d,%d) golij at j %d\n",
                    j,(unsigned int) (sums3D[k]>>(j1<<2))&0xf,(unsigned int) (sums3D[k]>>(j&0x3c))&0xf,(unsigned int) (sumsdebug3[k]>>(j1<<2))&0xf,
                    (unsigned int) (sums16[k]>>(j1<<2))&0xf,(unsigned int) (sumsdebug0[k]>>(j1<<2))&0xf,(unsigned int) (sumsdebug1[k]>>(j1<<2))&0xf,(unsigned int) (sumsdebug2[k]>>(j1<<2))&0xf,
                    sum1,(unsigned int) (sums16[j&0x3]>>(j&0x3c))&0xf,sum2,sum3,(unsigned int) (golij>>j)&1);
            }
        }
        //******************************************************************************
