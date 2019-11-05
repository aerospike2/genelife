//
//  external_set.c
//  project genelife
//
//  external setting of various quantities (typically called from python via genelife_c_interface.py)
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
//-------------------------------------------------------------------- set ...--------------------------------------------------------------------------
void set_colorfunction(int colorfunctionin) {
    if((colorfunctionin>12) || (colorfunctionin<0)) fprintf(stderr,"error colorfunction value passed %d out of range\n",colorfunctionin);
    else     colorfunction = colorfunctionin;
}
//.......................................................................................................................................................
void set_colorfunction2(int colorfunctionin) {
    if((colorfunctionin>12) || (colorfunctionin<-1)) fprintf(stderr,"error colorfunction value passed %d out of range\n",colorfunctionin);
    else     colorfunction2 = colorfunctionin;
}
//.......................................................................................................................................................
int setget_act_ymax(int actymax) {                  // sets ymax for activities only if argument nonzero, reads old value
    int ymaxold;
    ymaxold = ymax;
    ymax = actymax;
    return(ymaxold);
}
//.......................................................................................................................................................
int setget_act_ymaxq(int actymaxq) {                  // sets ymax for activities only if argument nonzero, reads old value
    int ymaxqold;
    ymaxqold = ymaxq;
    ymaxq = actymaxq;
    return(ymaxqold);
}
//.......................................................................................................................................................
void set_selectedgene(uint64_t gene) {
    selectedgene=gene;
    fprintf(stderr,"selected gene set to %llx\n",selectedgene);
}
//.......................................................................................................................................................
void set_offsets(int dx,int dy,int dt) {
    offdx =dx;
    offdy = dy;
    if(dt>0) {
        dt=0;
        fprintf(stderr,"positive time offsets not allowed, looking into the future not possible\n");
    }
    if(dt<=-maxPlane) {
        dt=-maxPlane+1;
        fprintf(stderr,"not enough planes set for this time offset, recompile software with larger maxPlane value\n");
    }
    offdt = dt;
}
//.......................................................................................................................................................
void set_quadrant(int quadrant) {
    if (quadrant >= -1 && quadrant < 7) quadrants = quadrant;
    repscheme &= ~R_quadrant;                                           // remove all quadrant bits : only one set at a time in interactive version
    if(quadrant >= 0 && quadrant < 7) {
        repscheme |= R_16_quadrant_sele<<quadrant;                          // assumes quadrant selectors are 7 successive bits following R_16_...
    }
}
//.......................................................................................................................................................
void set_randominflux(int randominfluxin) {
    randominflux=randominfluxin;
}
//.......................................................................................................................................................
void set_rbackground(int rbackgroundin, int randominfluxin) {
    rbackground=rbackgroundin;
    randominflux=randominfluxin;
    if(!rbackground) randominflux=0;   // do not leave patch randominflux variable active when turning off random background
}
//.......................................................................................................................................................
unsigned int set_repscheme_bits(int quadrant, int x, int y, unsigned int surviveover[]) {
    unsigned int quadrantval;

    quadrantval=(x<(Nmask>>1)? 0 : 1) + (y<(Nmask>>1)? 0 : 2);              // determine selected quadrant
    if(quadrants >= 0 && quadrants < 5) {                                   // assumes repscheme bits in pairs starting from bit 0 matching quadrants
        repscheme &=  ~(0x3llu<<(quadrants<<1));
        repscheme |=  quadrantval<<(quadrant<<1);
    }
    else if (quadrant < 6) {
        survivalmask &= ~0x3u;
        survivalmask|= quadrantval;
    }
    else if (quadrant<7) {
        overwritemask &= ~0x3u;
        overwritemask|= quadrantval;
    }
    repscheme &= ~R_quadrant;                                                // remove all quadrant bits : only one set at a time in interactive version
    quadrants = -1;                                                          // reset internal quadrants choice so that full display is shown
    surviveover[0]=survivalmask;
    surviveover[1]=overwritemask;
    return(repscheme);
}
//.......................................................................................................................................................
void set_repscheme(unsigned int repscheme_in) {
    repscheme = repscheme_in;
}
//.......................................................................................................................................................
void set_rulemod(unsigned int rulemod_in) {
    rulemod = rulemod_in;
}
//.......................................................................................................................................................
void set_surviveover64(unsigned int surviveover[], int len ) {
    if (len==3) {
        survivalmask = surviveover[0];
        if(selection<8) overwritemask = surviveover[1];
        else {
            birthmask = surviveover[1];
            overwritemask = surviveover[2];
        }
    }
    else fprintf(stderr,"surviveover64 needs three parameters, %d provided\n",len);
}
//.......................................................................................................................................................
void set_vscrolling() {
    vscrolling=1-vscrolling;
    if(!vscrolling) last_scrolled = 0;
}
//.......................................................................................................................................................
void set_noveltyfilter() {
    noveltyfilter=1-noveltyfilter;
}
//.......................................................................................................................................................
void set_activity_size_colormode() {
    activity_size_colormode = (activity_size_colormode+1) % 4;
}
//.......................................................................................................................................................
void set_gcolors() {
    gcolors = (gcolors+1)%10;
}
//.......................................................................................................................................................
void set_seed(int seed) {
    ranseed = seed;
}
//.......................................................................................................................................................
void set_nbhist(int nbhistin) {
    if(nbhist<nNhist*2) nbhist=nbhistin;
    else fprintf(stderr,"nbhist out of range %d > %d\n",nbhistin,nNhist*2-1);
}
//.......................................................................................................................................................
void set_genealogycoldepth(int genealogycoldepthin) {
    genealogycoldepth = genealogycoldepthin;
}
//.......................................................................................................................................................
void set_ancestortype(int ancestortypein) {
    if(ancestortypein <3) ancestortype = ancestortypein;
    else fprintf(stderr,"ancestor type %d out of range [0..2]\n",ancestortypein);
}
//.......................................................................................................................................................
void set_info_transfer_h(int do_info_transfer, int nbhood) {
    info_transfer_h = do_info_transfer;
    if(nbhood == 3 || nbhood == 5 || nbhood == 7)
        it_nbhood = nbhood;
    else fprintf(stderr,"error in nbhood value %d, allowed values are 3,5,7\n",nbhood);
}
//.......................................................................................................................................................
void set_activityfnlut(int activityfnlutin) {
    activityfnlut = activityfnlutin;
}
//.......................................................................................................................................................
void set_colorupdate1(int update1) {
    colorupdate1 = update1;
}
