//
//  maxmatch.c
//  fastgenegol
//
//  Created by John McCaskill on 21/01/2019.
//  Copyright © 2019 ECLT. All rights reserved.
//  https://www.geeksforgeeks.org/hopcroft-karp-algorithm-for-maximum-matching-set-2-implementation/
//  C++ implementation of Hopcroft Karp algorithm for maximum matching
//  adapted to handle sparse matrix as in LAPMOD

#include <stdio.h>

#define NIL 0
#define INF 0xffff   /* INT_MAX 0x7fffffff*/
#define TRUE 1
#define FALSE 0
typedef char boolean;

// C routines to implement one queue using an array
// if a number of new queues must be created, need constructor and destructor code via malloc, free
extern const int NLM;
int qrear = - 1;
int qfront = - 1;

extern inline void queue_push(int value)
{
    extern int queue_array[];
    if ((qfront == 0 && qrear == NLM-1) ||
            (qrear == (qfront-1)%(NLM-1)))
    {
        printf("Queue Overflow\n");
        return;
    }
  
    else if (qfront == -1) /* Insert First Element */
    {
        qfront = qrear = 0;
        queue_array[qrear] = value;
    }
  
    else if (qrear == NLM-1 && qfront != 0)
    {
        qrear = 0;
        queue_array[qrear] = value;
    }
  
    else
    {
        qrear++;
        queue_array[qrear] = value;
    }
}

extern inline int queue_pop()
{
    extern int queue_array[];
    if (qfront == - 1 || qfront > qrear) {
        fprintf(stderr,"Queue Underflow \n");
        return -1;
    }
    int data = queue_array[qfront];
    queue_array[qfront] = -1;
    if (qfront == qrear) {
        qfront = -1;
        qrear = -1;
    }
    else if (qfront == NLM-1)
        qfront = 0;
    else {
        // fprintf(stderr,"Element deleted from queue is : %d\n", data);
        qfront++;
    }
    return data;
}

// The graph is represented in sparse matrix form using cclap,kklap and iilap
// m and n are left and right hand side lists of vertices, NB 0 is used for dummy vertex here CHECK
// adj[u] storing adjacents of left side vertex 'u' is found from kk and ii
// Returns true if there is an augmenting path, else returns
// false
boolean bfs(int m, short unsigned int kk[], unsigned int ii[], short unsigned int pairU[], short unsigned int pairV[], short unsigned int dist[])
{
    // queue<int> Q; //an integer queue
  
    // First layer of vertices (set distance as 0)
    for (int u=1; u<=m; u++)
    {
        // If this is a free vertex, add it to queue
        if (pairU[u]==NIL)
        {
            // u is not matched
            dist[u] = 0;
            queue_push(u);
        }
  
        // Else set distance as infinite so that this vertex
        // is considered next time
        else dist[u] = INF;
    }
  
    // Initialize distance to NIL as infinite
    dist[NIL] = INF;
  
    // Q is going to contain vertices of left side only.
    while (qfront != - 1 && qfront <= qrear)
    {
        // Dequeue a vertex
        int u = queue_pop();
  
        // If this node is not NIL and can provide a shorter path to NIL
        if (dist[u] < dist[NIL])
        {
            // Get all adjacent vertices of the dequeued vertex u
            // list<int>::iterator i;
            // for (i=adj[u].begin(); i!=adj[u].end(); ++i)
            for (int i=ii[u-1]; i<ii[u]; i++)
            {
                int v = kk[i];
  
                // If pair of v is not considered so far
                // (v, pairV[V]) is not yet explored edge.
                if (dist[pairV[v]] == INF)
                {
                    // Consider the pair and add it to queue
                    dist[pairV[v]] = dist[u] + 1;
                    queue_push(pairV[v]);
                }
            }
        }
    }
  
    // If we could come back to NIL using alternating path of distinct
    // vertices then there is an augmenting path
    return (dist[NIL] != INF);
}

// Returns true if there is an augmenting path beginning with free vertex u
boolean dfs(int u, short unsigned int kk[], unsigned int ii[], short unsigned int pairU[], short unsigned int pairV[], short unsigned int dist[])
{
    if (u != NIL)
    {
        // list<int>::iterator i;
        // for (i=adj[u].begin(); i!=adj[u].end(); ++i)
        for (int i=ii[u-1]; i<ii[u]; i++)
        {
            // Adjacent to u
            // int v = *i;
            int v = kk[i];
  
            // Follow the distances set by BFS
            if (dist[pairV[v]] == dist[u]+1)
            {
                // If dfs for pair of v also returns
                // true
                if (dfs(pairV[v],kk,ii,pairU,pairV,dist) == TRUE)
                {
                    pairV[v] = u;
                    pairU[u] = v;
                    return TRUE;
                }
            }
        }
  
        // If there is no augmenting path beginning with u.
        dist[u] = INF;
        return FALSE;
    }
    return TRUE;
}

int hopcroftKarp(int m, int n, short unsigned int kk[], unsigned int ii[], short unsigned int pairU[], short unsigned int pairV[], short unsigned int dist[]) { // Returns size of maximum matching
    // pairU[u] stores pair of u in matching where u is a vertex on left side of Bipartite Graph.
    // If u doesn't have any pair, then pairU[u] is NIL
    // pairV[v] stores pair of v in matching. If v doesn't have any pair, then pairU[v] is NIL
    // dist[u] stores distance of left side vertices
    // dist[u] is one more than dist[u'] if u is next to u'in augmenting path
  
    for (int u=0; u<=m; u++)                     // Initialize NIL as pair of all vertices
        pairU[u] = NIL;
    for (int v=0; v<=n; v++)                     // Note that not all of right vertices n have edges to them
        pairV[v] = NIL;
  
    int result = 0;                              // Initialize result
  
    while (bfs(m,kk,ii,pairU,pairV,dist))        // Keep updating the result while there is an augmenting path.
    {
        for (int u=1; u<=m; u++)                 // Find a free vertex
                                                 // If current vertex is free and there is an augmenting path from current vertex
            if (pairU[u]==NIL && dfs(u,kk,ii,pairU,pairV,dist))
                result++;
    }
    return result;
}

int maxmatch(int m, short unsigned int kk[], unsigned int ii[], short unsigned int pairU[], short unsigned int pairV[], short unsigned int dist[])
{
    int i,n;
    for (n=i=0;i<ii[m];i++) n = kk[i] > n ? kk[i] : n;    // calculate n as maximum vertex on right
    i=hopcroftKarp(m,n,kk,ii,pairU,pairV,dist);
    // fprintf(stderr,"Size of maximum matching is %d\n",i);
    return i;
}
