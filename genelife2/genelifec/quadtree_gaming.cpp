///////////////////////////////////////////////////////////////////////////
//  This file is part of Quicksilver - a bike messenger simulation game  //
//  Copyright (C) 2005 Scott Czepiel <http://czep.net/>                  //
//                                                                       //
//  This program is free software; you can redistribute it and/or modify //
//  it under the terms of the GNU General Public License as published by //
//  the Free Software Foundation; either version 2 of the License, or    //
//  (at your option) any later version.                                  //
//                                                                       //
//  This program is distributed in the hope that it will be useful,      //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of       //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        //
//  GNU General Public License for more details.                         //
//                                                                       //
//  You should have received a copy of the GNU General Public License    //
//  along with this program; if not, write to the Free Software          //
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 //
//  USA                                                                  //
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
//                                                                       //
// quadtree.hpp - header file for PMR Quadtree class                     //
//                                                                       //
// References:                                                           //
//                                                                       //
// Erik Hoel and Hanan Samet, "Efficient Processing of Spatial           //
//   Queries in Line Segment Databases", Advances in Spatial             //
//   Databases - 2nd Symp., SSD '91, (O. Gunther and H. J. Schek,        //
//   eds.), Lecture Notes in Computer Science 525, Springer-Verlag,      //
//   Berlin, 1991, 237-256.                                              //
//                                                                       //
// G. R. Hjaltason and H. Samet, "Ranking in spatial databases", in      //
//   Advances in Spatial Databases - 4th Symposium, SSD'95, M. J.        //
//   Egenhofer and J. R. Herring, Eds., Lecture Notes in Computer        //
//   Science 951, Springer-Verlag, Berlin, 1995, 83-95.                  //
//                                                                       //
// Gisli Hjaltason and Hanan Samet, "Speeding Up Construction of PMR     //
//   Quadtree-Based Spatial Indexes," The VLDB Journal, 11 (2002) 2,     //
//   pp. 109-137. Springer-Verlag.                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////


#include "quadtree.hpp"


// MAX and MIN macros used in EdgeInteresects and BoxIntersects functions
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))


//
// class of elements in the NearestEdge priority queue
//
class PQElem 
{
friend class PMRQuadtree;
private:
    bool IsEdge;    // true if elem represents an edge, false if a PMRBlock
    double dist;    // distance to query point, the key of the PQ
    PMRBlock* Block;// pointer to PMRBlock if IsEdge==false
    int edge_idx;   // index in E of edge if IsEdge==true
public:
    PQElem(bool isedge = true, double d = -1.0, PMRBlock* b = 0, int eid = -1)
        : IsEdge(isedge), dist(d), Block(b), edge_idx(eid) {}
    friend bool operator<(const PQElem& x, const PQElem& y)
    {
        return (x.dist > y.dist);
    }
    void operator=(const PQElem& RHS)
    {
        IsEdge = RHS.IsEdge;
        dist = RHS.dist;
        Block = RHS.Block;
        edge_idx = RHS.edge_idx;
    }
};        


//
// PMRBlock constructor
//
PMRBlock::PMRBlock(int set_depth)
{
    _depth = set_depth;
    _parent = 0;
    _NE = 0;
    _NW = 0;
    _SW = 0;
    _SE = 0;
    _x = 0;
    _y = 0;
}


//
// PMRQuadtree constructor
//
PMRQuadtree::PMRQuadtree(int set_depth, int set_k)
{
    _depth = set_depth;
    _k = set_k;
    _root = 0;
    _numleaves = 0;
}

//
// Build - build PMR Quadtree from the Graph g
//
int PMRQuadtree::Build()
{   
    int i, x1, y1, x2, y2;
    PMRBlock* Block;
    std::queue<PMRBlock*> Q;
    
    // assert that we have a valid pointer to graph
    if (_g == 0) return 1;
    
    // setup the root block of the quadtree
    try {
        _root = new PMRBlock(_depth);
    }
    catch (std::bad_alloc) {
        return 1;
    }

    //
    // load each segment into the tree
    //    
    for (i = 0; i < _g->Ecount; i++) {
        
        // assert that the queue is initially empty
        if (Q.size() != 0) {
            return 1;
        }
    
        // get the coordinates of the current segment
        x1 = _g->E.x[0][i];
        y1 = _g->E.y[0][i];
        x2 = _g->E.x[1][i];
        y2 = _g->E.y[1][i];
    
        // enqueue the root of the tree
        Q.push(_root);
        
        // do while blocks are still queued
        while (Q.size() > 0) {
            
            // pop the next block from the queue
            Block = Q.front();
            
            // do if this is a leaf node
            if (Block->_NE == 0) {                
                
                // test if edge intersects with or is contained within block
                if (EdgeIntersects(Block, x1, y1, x2, y2)) {
                    
                    // add edge to the list for this block
                    Block->_e.push_back(i);
                    
                    // split the block if we exceed the threshold
                    if (Block->_e.size() > _k) {
                        if (Split(Block) != 0) {
                            return 1;
                        }
                    }
                }                
            }
                        
            // do if non-leaf node
            else {
                
                // test each of the four children against the bounding
                // box of the current edge, enqueue any matching blocks
                
                // NE
                if (BoxIntersects(Block->_NE, x1, y1, x2, y2)) {
                    Q.push(Block->_NE);
                }
                // NW
                if (BoxIntersects(Block->_NW, x1, y1, x2, y2)) {
                    Q.push(Block->_NW);
                }
                // SW
                if (BoxIntersects(Block->_SW, x1, y1, x2, y2)) {
                    Q.push(Block->_SW);
                }
                // SE
                if (BoxIntersects(Block->_SE, x1, y1, x2, y2)) {
                    Q.push(Block->_SE);
                }
                
            }
            // end switch on leaf or non-leaf node
            
            // dequeue the current block
            Q.pop();
            
        }
        // end while blocks still exist in Q        
    }
    // end of loop for each segment
       
    return 0;
} 


//
// Split - split a block into 4 child blocks
//    return 0 on success, non-zero on something bad
int PMRQuadtree::Split(PMRBlock* Block)
{   
    int         i, j, ii;
    int         x1, y1, x2, y2;
    PMRBlock*   NewBlock;
            
    // fail if block already has children
    if (Block->_NE != 0)
        return 1;
        
    // return silently if depth has already been reached
    if (Block->_depth == 0)
        return 0;
    
    // allocate space for new child blocks
    try {
        Block->_NE = new PMRBlock(Block->_depth - 1);
        Block->_NW = new PMRBlock(Block->_depth - 1);
        Block->_SW = new PMRBlock(Block->_depth - 1);
        Block->_SE = new PMRBlock(Block->_depth - 1);
    }
    catch (std::bad_alloc) {
        return 1;
    }

    // set each new child block's parent
    Block->_NE->_parent = Block;
    Block->_NW->_parent = Block;
    Block->_SW->_parent = Block;
    Block->_SE->_parent = Block;
    
    // set each child block's origin coordinates 
    Block->_NE->_x = Block->_x + (2 << Block->_NE->_depth);
    Block->_NE->_y = Block->_y + (2 << Block->_NE->_depth);
    
    Block->_NW->_x = Block->_x;
    Block->_NW->_y = Block->_y + (2 << Block->_NE->_depth);
    
    Block->_SW->_x = Block->_x;
    Block->_SW->_y = Block->_y;
    
    Block->_SE->_x = Block->_x + (2 << Block->_NE->_depth);
    Block->_SE->_y = Block->_y;

    // increment the number of tree leaves by 3,
    // remember, previously the parent was considered a leaf 
    _numleaves += 3;
    
    //
    // Build the child edge lists from the edges in the parent list
    //
    
    for (i = 0; i < Block->_e.size(); i++) {
        
        // get the index in E of this edge in e
        ii = Block->_e[i];
        
        // start with the NE child block
        for (j = 0, NewBlock = Block->_NE; j < 4; j++) {
                                    
            x1 = _g->E.x[0][ii];
            y1 = _g->E.y[0][ii];
            x2 = _g->E.x[1][ii];
            y2 = _g->E.y[1][ii];
            
            // test if segment intersects or is entirely contained within block
            if (EdgeIntersects(NewBlock, x1, y1, x2, y2)) {
                
                // add edge to the list for new block
                NewBlock->_e.push_back(ii);  
                
                // NOTE: do not split the child block here even if it
                // contains all the elements of its parent, only split
                // the next time a segment is added to it

            }
            
            // then go to NW, SW, and SE blocks
            if (i == 0)
                NewBlock = Block->_NW;
            else if (i == 1)
                NewBlock = Block->_SW;
            else if (i == 2)
                NewBlock = Block->_SE;
                
        }   // continue to next child block
  
    }   // continue to next edge in parent's original list
    
    // now remove all items from the parent block's list
    Block->_e.clear();
        
    return 0;
    
}

//
// EdgeIntersects - determine whether a given segment intersects with 
//   a block, or alternatively if the segment is entirely contained 
//   within the block
//
int PMRQuadtree::EdgeIntersects(PMRBlock* Block, int x1, int y1, int x2, int y2)
{   
    // This algorithm adapted from Kyle Louden's "Mastering Algorithms in C" 

    int     x3, x4, y3, y4;
    int     z1, z2, s1, s2;
    int     i, w, match;
    
    // define the block as a bounding box 
    w = (2 << Block->_depth);
    x3 = Block->_x;
    y3 = Block->_y;
    x4 = Block->_x + w;
    y4 = Block->_y + w;
    
    // first test whether the segment bounding box is entirely contained
    // within the block.  If true, the segment is obviously inside the block
    
    if ((x3 <= MIN(x1, x2)) &&
        (MAX(x1, x2) <= x4) &&
        (y3 <= MIN(y1, y2)) &&
        (MAX(y1, y2) <= y4)) return 1;
        
    // if segment bounding box is not entirely contained within the block, then
    // the segment itself must intersect with at least one of the four sides
    // of the block.  Otherwise, the segment and block do not intersect 
    
    // reset y4 such that we are testing the south edge of the block 
    y4 = Block->_y;
    
    for (i = 0, match = 0; i < 4; i++) {
    
        // determine orientation of p3 relative to p2, wrt p1 
        if ((z1 = ((x3-x1)*(y2-y1)) - ((y3-y1)*(x2-x1))) < 0)
            s1 = -1;
        else if (z1 > 0)
            s1 = 1;
        else
            s1 = 0;
        
        // determine orientation of p4 relative to p2, wrt p1 
        if ((z2 = ((x4-x1)*(y2-y1)) - ((y4-y1)*(x2-x1))) < 0)
            s2 = -1;
        else if (z2 > 0)
            s2 = 1;
        else
            s2 = 0;
    
        // if the signs of z1 and z2 are different, or either is 0, 
        // then signal a match 
        if ((s1 == 0 || s2 == 0) || (s1 != s2)) {
            match = 1;
            break;
        }
        
        // prepare to test another edge
        
        // 2nd test = west edge 
        if (i == 0) {
            x4 = Block->_x;
            y4 = Block->_y + w;
        }
        
        /* 3rd test = north edge */
        else if (i == 1) {
            x3 = Block->_x + w;
            y3 = Block->_y + w;
        }
        
        /* 4th test = east edge */
        else if (i == 2) {
            x4 = Block->_x + w;
            y4 = Block->_y;
        }
        
    }
    
    if (match)
        return 1;
    else    
        return 0;
        
}


//
// BoxIntersects - determine whether a segment bounding box 
//   intersects with block
//
int PMRQuadtree::BoxIntersects(PMRBlock* Block, int x1, int y1, int x2, int y2)
{
    
    int x3, x4, y3, y4;
    int w;
        
    // define the block as a bounding box 
    w = (2 << Block->_depth);
    x3 = Block->_x;
    y3 = Block->_y;
    x4 = Block->_x + w;
    y4 = Block->_y + w;
    
    // perform the test, return 1 if the boxes intersect 
    if ((MAX(x1, x2) >= MIN(x3, x4)) &&
        (MAX(x3, x4) >= MIN(x1, x2)) &&
        (MAX(y1, y2) >= MIN(y3, y4)) &&
        (MAX(y3, y4) >= MIN(y1, y2))) return 1;
    

    return 0;        
}


//
// NearestEdge - find the segment nearest to a given query point
//   return the index in G->E of the nearest segment if found,
//   else return -1 
int PMRQuadtree::NearestEdge(int qx, int qy, double radius)
{
    int retval = -1;
    bool found = false;
    double dist;
    int i;
    int edge_count;
    int eid;
    int x1, x2, y1, y2;
    PMRBlock* Block;
    PQElem pqelem;
    std::priority_queue<PQElem> PQ;

    // quick return if query point is out of range
    if (qx < 0 || qx >= (2 << _depth) || qy < 0 || qy >= (2 << _depth))
        return retval;
    
    // if radius parameter is 0, set to infinite
    if (radius == 0.0) radius = DBL_MAX;
    
    // set the root block of the quadtree as the first
    // element of the priority queue
    PQ.push(PQElem(false, 0.0, _root, -1));
    
    //
    // main loop
    //
    while (!found && !PQ.empty()) {
        
        // pop the top of the priority queue, this element is either
        // an edge or a Block, in either case it is the nearest element
        // to the query point
        pqelem = PQ.top();
        PQ.pop();
        
        // if this is a segment, return its index in E
        // unless it is beyond the search radius
        if (pqelem.IsEdge) {
            
            found = true;
            if (pqelem.dist > radius) {
                retval = -1;
            }
            else {
                retval = pqelem.edge_idx;
            }
        }
        
        // if leaf block, then add its segments to the pq
        else if (pqelem.Block->_NE == 0) {
            
            edge_count = pqelem.Block->_e.size();
            
            for (i = 0; i < edge_count; i++) {
                
                eid = pqelem.Block->_e[i];
                x1 = _g->E.x[0][eid];
                y1 = _g->E.y[0][eid];
                x2 = _g->E.x[1][eid];
                y2 = _g->E.y[1][eid];
                
                // calculate dsq from query point to this segment
                dist = DistPointSegment(qx, qy, x1, y1, x2, y2);
                
                // add to PQ if not less than distance to current block
                // and less than radius
                if (dist >= pqelem.dist && dist <= radius) {                    
                    PQ.push(PQElem(true, dist, 0, eid));
                }
            }
        }
        
        // if non-leaf block, add its 4 children to the pq
        else {
            
            // start with the NE child
            Block = pqelem.Block->_NE;
            
            // loop 4 times
            for (i = 0; i < 4; i++) {
                
                // calculate the distance to this block
                dist = DistPointBlock(qx, qy, Block);
                
                // add the block to the pq if dist is less than radius
                if (dist <= radius) {
                    PQ.push(PQElem(false, dist, Block, -1));
                }
                
                // set to the next child block
                if (i == 0)
                    Block = pqelem.Block->_NW;
                else if (i == 1)
                    Block = pqelem.Block->_SW;
                else if (i == 2)
                    Block = pqelem.Block->_SE;
                    
            }
            
        }
        // END OF IF STATEMENT BASED ON TYPE OF ELEM IN THE PQ
        
    }
    // END OF MAIN LOOP WHILE NOT FOUND
    
    // now free any remaining elems in the pq
    // to be nice to the memory manager
    while (!PQ.empty()) {
        PQ.pop();
    }
        
    return retval;
}


//
// DistPointSegment - calculate the squared distance between a point and a line segment
//
double PMRQuadtree::DistPointSegment(int qx, int qy, int segx1, int segy1, int segx2, int segy2)
{    
    // for details, see the comp.graphics.algorithms FAQ
    
    double  x, y, Ax, Ay, Bx, By;
    double  d, r, s;
    
    x = (double) qx;
    y = (double) qy;
    
    Ax = (double) segx1;
    Bx = (double) segx2;
    Ay = (double) segy1;
    By = (double) segy2;
    
    // calculate squared length of the segment 
    d = (Bx-Ax)*(Bx-Ax) + (By-Ay)*(By-Ay);
    
    // if distance is zero, then the segment is just a point
    if (d == 0.0) {
        return ((x-Ax)*(x-Ax) + (y-Ay)*(y-Ay));
    }
    
    // calculate the parameter, r
    r = ((x-Ax)*(Bx-Ax) + (y-Ay)*(By-Ay)) / d;
    
    // r helps find the perpendicular from segment to query point
    // r < 0 : P falls to the left of segment, return distance to A
    // r = 0 : P = A, return distance to A
    // 0<r<1 : P is on the segment, return length of perp. segment
    // r = 1 : P = B, return distance to B
    // r > 1 : P falls to the right of segment, return distance to B 
    
    if (r <= 0) {
        return ((x-Ax)*(x-Ax) + (y-Ay)*(y-Ay));
    }
    else if (r >= 1) {
        return ((x-Bx)*(x-Bx) + (y-By)*(y-By));
    }
    else {
        s = ((Ay-y)*(Bx-Ax)-(Ax-x)*(By-Ay));
        return (s * s / d);
    }    
}


//
// DistPointBlock - calculate the squared distance between a point and a PMRBlock
// note: if block contains point, return 0
double PMRQuadtree::DistPointBlock(int x, int y, PMRBlock* Block)
{
    int     w;
    double  dist = 0.0;
    double  dnw, dne, dsw, dse, min, min2;
    int     v1[3], v2[3], v3[3], v4[3];
    
    // store the width of the block
    w = (2 << Block->_depth);
    
    // if block contains point, return 0
    if (Block->_x <= x &&
        x < Block->_x + w &&
        Block->_y <= y &&
        y < Block->_y + w) return dist;
    
    // calculate the distance from point to each of the corners of the block 
    v1[1] = Block->_x;
    v1[2] = Block->_y;
    v2[1] = Block->_x;
    v2[2] = Block->_y;
    v3[1] = Block->_x;
    v3[2] = Block->_y;
    v4[1] = Block->_x;
    v4[2] = Block->_y;
    
    
    dsw = (x-Block->_x)*(x-Block->_x) + (y-Block->_y)*(y-Block->_y);
    dnw = (x-Block->_x)*(x-Block->_x) + (y-(Block->_y+w))*(y-(Block->_y+w));
    dse = (x-(Block->_x+w))*(x-(Block->_x+w)) + (y-Block->_y)*(y-Block->_y);
    dne = (x-(Block->_x+w))*(x-(Block->_x+w)) + (y-(Block->_y+w))*(y-(Block->_y+w));
    
    // find the nearest point 
    v1[1] = Block->_x;
    v1[2] = Block->_y;
    min = dsw;
    if (dnw < min) {
        v1[1] = Block->_x;
        v1[2] = Block->_y + w;
        min = dnw;
    }
    if (dse < min) {
        v1[1] = Block->_x + w;
        v1[2] = Block->_y;
        min = dse;
    }
    if (dne < min) {
        v1[1] = Block->_x + w;
        v1[2] = Block->_y + w;
        min = dne;
    }
    
    // find the second nearest point
    v2[1] = Block->_x;
    v2[2] = Block->_y;
    min2 = dsw;
    if (min2 == min) {
        v2[1] = Block->_x;
        v2[2] = Block->_y + w;
        min2 = dnw;
    }
    if (dse < min2 && dse != min) {
        v2[1] = Block->_x + w;
        v2[2] = Block->_y;
        min2 = dse;
    }
    if (dne < min2 && dne != min) {
        v2[1] = Block->_x + w;
        v2[2] = Block->_y + w;
        min2 = dne;
    }
    
    // get the distance to the segment formed by the two nearest points
    dist = DistPointSegment(x, y, v1[1], v1[2], v2[1], v2[2]);
    
    return dist;
}

