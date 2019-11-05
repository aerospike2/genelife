//
//  stash.c
//  project genelife
//
//  miscellaneous functions associated with stashing, saving and retrieving data
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
//---------------------------------------------------------------- save and retrieve gol... data --------------------------------------------------------
//.......................................................................................................................................................
void stash(){               // stash current gol,golg
    int ij;
    for (ij=0; ij<N2; ij++) {
        stashgol[ij] = planes[curPlane][ij];
        stashgolg[ij] = planesg[curPlane][ij];
        stashgolb[ij] = planesb[curPlane][ij];
        stashgolr[ij] = planesr[curPlane][ij];
        stashgolgstats[ij] = planesgs[curPlane][ij];
    }
}
//.......................................................................................................................................................
void label2stash(int cumul) { // stash current gol,golg, golb, golr, golgstats from selected labelled component (either cumulatively or individually)
    int ij;
    for (ij=0; ij<N2; ij++) {
        if (labelcc[ij]==0xffffffff) {
            stashgol[ij] = planes[curPlane][ij];
            stashgolg[ij] = planesg[curPlane][ij];
            stashgolb[ij] = planesb[curPlane][ij];
            stashgolr[ij] = planesr[curPlane][ij];
            stashgolgstats[ij] = planesgs[curPlane][ij];
        }
        else if (!cumul) {
            stashgol[ij] = 0ull;
            stashgolg[ij] = 0ull;
            stashgolb[ij] = 0ull;
            stashgolr[ij] = 0ull;
            stashgolgstats[ij] = 0ull;
        }
    }
}
//.......................................................................................................................................................
void unstash(){               // retrieve current gol,golg from stashed values
    int ij;
    for (ij=0; ij<N2; ij++) {
        planes[curPlane][ij] = stashgol[ij];
        planesg[curPlane][ij] = stashgolg[ij];
        planesb[curPlane][ij] = stashgolb[ij];
        planesr[curPlane][ij] = stashgolr[ij];
        planesgs[curPlane][ij] = stashgolgstats[ij];
    }
}
//.......................................................................................................................................................
int savegols( int step, uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[],uint64_t golr[]) {
    FILE *fp;
    char fname[30];

    sprintf(fname,"gol_gsbr_data%d.ext",step);
    fp = fopen( fname , "wb" );
    
    fwrite(gol  , sizeof(uint64_t) , N2 , fp );
    fwrite(golg , sizeof(uint64_t) , N2 , fp );
    fwrite(golgstats, sizeof(uint64_t) , N2 , fp );
    fwrite(golb , sizeof(uint64_t) , N2 , fp );
    fwrite(golr , sizeof(uint64_t) , N2 , fp );

    fclose(fp);
    return 0;
}
//---------------------------------------------------------------- save data ----------------------------------------------------------------------------
int retrievegols( int step, uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[],uint64_t golr[]) {
    FILE *fp;
    char fname[30];

    sprintf(fname,"gol_gsbr_data%d.ext",step);
    fp = fopen( fname , "rb" );
    
    fread(gol ,  sizeof(uint64_t) , N2, fp );
    fread(golg ,  sizeof(uint64_t) , N2, fp );
    fread(golgstats, sizeof(uint64_t) , N2, fp );
    fread(golb , sizeof(uint64_t) , N2, fp );
    fread(golr , sizeof(uint64_t) , N2, fp );

    fclose(fp);
    return 0;
}


