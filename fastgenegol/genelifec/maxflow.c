//
//  maxflow.c
//  genelife
//
// C++ program for implementation of Ford Fulkerson algorithm
// Adapted from https://www.geeksforgeeks.org/ford-fulkerson-algorithm-for-maximum-flow-problem/
// Simple implementation of queue

#include <stdio.h>
#include <limits.h>
#include <string.h>
#include "lapjv.h"

// Number of vertices in given graph
const int V=6;
const int QMAX=8*8;

/* Returns true if there is a path from source 's' to sink 't' in
residual graph. Also fills parent[] to store the path */

// C routines to implement one queue using an array
// if a number of new queues must be created, need constructor and destructor code via malloc, free

int qrear = - 1;
int qfront = - 1;
int queue_array[QMAX];

extern inline void queue_push(int value)
{
    if ((qfront == 0 && qrear == QMAX-1) ||
            (qrear == (qfront-1)%(QMAX-1)))
    {
        printf("Queue Overflow\n");
        return;
    }
  
    else if (qfront == -1) /* Insert First Element */
    {
        qfront = qrear = 0;
        queue_array[qrear] = value;
    }
  
    else if (qrear == QMAX-1 && qfront != 0)
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
    else if (qfront == QMAX-1)
        qfront = 0;
    else {
        // fprintf(stderr,"Element deleted from queue is : %d\n", data);
        qfront++;
    }
    return data;
}

boolean bfs(int rGraph[V][V], int s, int t, int parent[])
{
    // Create a visited array and mark all vertices as not visited
    boolean visited[V];
    
    memset(visited, 0, sizeof(visited));
    
    // Create a queue, enqueue source vertex and mark source vertex as visited
    queue_push(s);
    visited[s] = TRUE;
    parent[s] = -1;

    // Standard BFS Loop
    while (qfront != - 1 && qfront <= qrear)
    {
        int u = queue_pop();

        for (int v=0; v<V; v++)
        {
            if (visited[v]==FALSE && rGraph[u][v] > 0)
            {
                queue_push(v);
                parent[v] = u;
                visited[v] = TRUE;
            }
        }
    }

    // If we reached sink in BFS starting from source, then return
    // true, else false
    return (visited[t] == TRUE);
}

// Returns the maximum flow from s to t in the given graph
int fordFulkerson(int graph[V][V], int s, int t)
{
    int u, v;

    // Create a residual graph and fill the residual graph with
    // given capacities in the original graph as residual capacities
    // in residual graph
    int rGraph[V][V]; // Residual graph where rGraph[i][j] indicates
                    // residual capacity of edge from i to j (if there
                    // is an edge. If rGraph[i][j] is 0, then there is not)
    for (u = 0; u < V; u++)
        for (v = 0; v < V; v++)
            rGraph[u][v] = graph[u][v];

    int parent[V]; // This array is filled by BFS and to store path

    int max_flow = 0; // There is no flow initially

    // Augment the flow while tere is path from source to sink
    while (bfs(rGraph, s, t, parent))
    {
        // Find minimum residual capacity of the edges along the
        // path filled by BFS. Or we can say find the maximum flow
        // through the path found.
        int path_flow = INT_MAX;
        for (v=t; v!=s; v=parent[v])
        {
            u = parent[v];
            path_flow = path_flow < rGraph[u][v]? path_flow : rGraph[u][v];
        }

        // update residual capacities of the edges and reverse edges
        // along the path
        for (v=t; v != s; v=parent[v])
        {
            u = parent[v];
            rGraph[u][v] -= path_flow;
            rGraph[v][u] += path_flow;
        }

        // Add path flow to overall flow
        max_flow += path_flow;
    }

    // Return the overall flow
    return max_flow;
}

// Driver program to test above functions
void testflow()
{
    // Let us create a graph shown in the above example
    int graph[V][V] = { {0, 16, 13, 0, 0, 0},
                        {0, 0, 10, 12, 0, 0},
                        {0, 4, 0, 0, 14, 0},
                        {0, 0, 9, 0, 0, 20},
                        {0, 0, 0, 7, 0, 4},
                        {0, 0, 0, 0, 0, 0}
                    };
    fprintf(stderr,"The maximum possible flow is %d\n",fordFulkerson(graph, 0, 5));

}

