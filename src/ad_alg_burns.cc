//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
//
// An implementation of Burns's minimum mean cycle algorithm.
//
// @techreport{Bu91,
// author = {S.~M. Burns},
// title = {Performance Analysis and Optimization of Asynchronous Circuits},
// institution = {California Institute of Technology},
// year = 1991,
// type = {{PhD} Thesis},
// address = "",
// number = "",
// month = "",
// }
//

#include "ad_graph.h"
#include "ad_queue.h"

// More node info for Burns's algorithm.
struct ninfo_burns {
    float dist;    // node distance or potential.
    int   length;  // path length from source in topological order.
    int   indeg;   // indegree.
};

#if 0
bool search( const ad_graph< ninfo > *g, int u, bool *visited, bool *critical )
{
    visited[ u ] = true;
    for ( int i = 0; i < g->outdegree( u ); ++i ) {
        int e = g->ith_target_edge( u, i );
        if ( critical[ e ] ) {
            int v = g->ith_target_node( u, i );
            if ( visited[ v ] )
                return true;
            else
                return search( g, v, visited, critical );
        }
    }
    return false;
}
#endif

float
find_min_cycle_mean_for_scc( const ad_graph< ninfo > *g, 
                             int plus_infinity,
                             float lambda_so_far )
{
    int n = g->num_nodes();
    int m = g->num_edges();

    // critical[e] is true if edge e is critical.

    ninfo_burns      *more_ninfo = new ninfo_burns[ n ];
    bool             *critical = new bool[ m ];
    ad_queue< int >  nodeq( n );

    float f_plus_infinity = ( float ) plus_infinity;
    float f_minus_infinity = -f_plus_infinity;

    // STEP: Initialize lambda to the minimum of the min edge weight and
    // the previous lambda:
    float lambda = lambda_so_far;
    for ( int e = 0; e < m; ++e )
        min2( lambda, ( float ) g->edge_info( e ) );

    // STEP: Initialize the distance of each node to zero:
    for ( int v = 0; v < n; ++v )
        more_ninfo[ v ].dist = 0.0;

    // STEP: Iterate until the critical graph is cyclic.
    while ( true ) {

#ifdef PROGRESS
        printf( "PROGRESS Iteration number= %d lambda= %10.2f\n", NITER, lambda );
#endif

#ifdef PROGRESS
        ++NITER;
#endif

        // STEP: Find all critical edges:
        int ncrit = 0;
        for ( int e = 0; e < m; ++e ) {
            int u = g->source( e );
            int v = g->target( e );

            float delta1 = more_ninfo[ u ].dist + g->edge_info( e ) - more_ninfo[ v ].dist;
            if ( fabs_val( lambda - delta1 ) < SMALL_EPSILON ) {
                critical[ e ] = true;
                ++ncrit;
            } else {
                critical[ e ] = false;
            }
        }  // for each edge

        // STEP: Topologically sort the critical graph:
    
        // Find the indegree and ( negative ) length of each node and put
        // them into the nodeq.
        nodeq.init();

        for ( int v = 0; v < n; ++v ) {
            more_ninfo[ v ].indeg = 0;

            for ( int i = 0; i < g->indegree( v ); ++i ) {
                if ( critical[ g->ith_source_edge( v, i ) ] )
                    more_ninfo[ v ].indeg++;
            }

            if ( 0 == more_ninfo[ v ].indeg ) {
                more_ninfo[ v ].length = 0;
                nodeq.put( v );
            } else {
                more_ninfo[ v ].length = plus_infinity;
            }
        }  // for

        // Do the actual topological sorting using the previously found
        // indegrees and lengths for nodes:

        int count_visited = 0;  // Number of visited nodes.
        while ( nodeq.is_not_empty() ) {
            ++count_visited;

            int u = nodeq.get();
            for ( int i = 0; i < g->outdegree( u ); ++i ) {
                if ( critical[ g->ith_target_edge( u, i ) ] ) {
                    // edge = u->v
                    int v = g->ith_target_node( u, i );

                    // length[ v ] = min( length[ v ], length[ u ] - 1 )
                    min2( more_ninfo[ v ].length, more_ninfo[ u ].length - 1 );
                    more_ninfo[ v ].indeg--;

                    if ( 0 == more_ninfo[ v ].indeg )
                        nodeq.put( v );
                }
            }
        }  // while

        // STEP: If the critical graph is cyclic, then the optimum lambda
        // is found, so exit.
        if ( count_visited != n )
            break;

#if 0
        /***********************/
        if ( count_visited != n ) {
            bool *visited = new bool [ n ];
            for ( int u = 0; u < n; ++u )
                visited[ u ] = false;
            bool found = false;
            for ( int u = 0; u < n; ++u ) {
                if ( ! visited[ u ] )
                    found = search( g, u, visited, critical );
            }
            if ( found )
                printf( "cycle found\n" );
            else
                printf( "no cycle found\n" );
            delete [] visited;
            break;
        }
        /***********************/
#endif

        // STEP: Find theta to update lambda as well as the distance of
        // every node:
        float theta = f_minus_infinity;
        for ( int e = 0; e < m; ++e ) {
            // e = u->v
            int u = g->source( e );
            int v = g->target( e );

            int delta2 = more_ninfo[ v ].length + 1 - more_ninfo[ u ].length;
            if ( delta2 > 0 ) {
                float delta1 = more_ninfo[ u ].dist + g->edge_info( e ) - more_ninfo[ v ].dist;
                max2( theta, ( lambda - delta1 ) / delta2 );
            }
        }  // for each edge

        // STEP: Using theta, update lambda as well as the distance of
        // every node:
        lambda -= theta;
        for ( int v = 0; v < n; ++v )
            more_ninfo[ v ].dist -= theta * more_ninfo[ v ].length;

    }  // main while loop

    delete [] more_ninfo;
    delete [] critical;

    return lambda;
} // find_min_cycle_mean_for_scc

// End of file
