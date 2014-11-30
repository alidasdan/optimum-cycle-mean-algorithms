//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
//
// An implementation of Lawler's minimum mean cycle algorithm using
// Bellman-Ford's shortest path algorithm.
//

#include "ad_graph.h"
#include "ad_cqueue.h"

// More node info for Lawler's algorithm.
struct ninfo_lawler {
    float dist;
    int   not_included;
};

float
find_min_cycle_mean_for_scc( const ad_graph< ninfo > *g, 
                             int plus_infinity,
                             float lambda_so_far )
{
    int n = g->num_nodes();
    int m = g->num_edges();

    float f_plus_infinity = ( float ) plus_infinity;
    float f_minus_infinity = -f_plus_infinity;

    // STEP: Determine lower and upper bounds on lambda. 
#ifdef IMPROVE_LAMBDA_BOUNDS
    float upper = find_min_lambda( g, plus_infinity );
    float lower = f_plus_infinity;

    for ( int e = 0; e < m; ++e ) {
        min2( lower, ( float ) g->edge_info( e ) );
    }
#else
    float upper = -f_plus_infinity;
    float lower = f_plus_infinity;

    for ( int e = 0; e < m; ++e ) {
        min2( lower, ( float ) g->edge_info( e ) );
        max2( upper, ( float ) g->edge_info( e ) );
    }
#endif

    if ( lambda_so_far <= lower )
        return lambda_so_far;

    min2( upper, ( float ) 2.0 * lambda_so_far - lower );

    float lambda = upper;

    ninfo_lawler      *more_ninfo = new ninfo_lawler[ n ];
    ad_cqueue< int >  nodeq( n + 1 );  // +1 for insertion of END_PHASE node.

#define END_PHASE -1

    while ( ( upper - lower ) > EPSILON ) {

        // Determine the new lambda in the middle.
        lambda = ( upper + lower ) / 2;

#ifdef PROGRESS
        ++NITER;
#endif

#ifdef PROGRESS
        printf( "PROGRESS Iteration number= %d lambda= %10.2f\n", NITER, lambda );
#endif

        // Subtract lambda from each edge weight.
        // Check to see if the resulting graph has a negative cycle.

        more_ninfo[ SOURCE ].dist = 0;
        more_ninfo[ SOURCE ].not_included = 0;
        for ( int v = 1; v < n; ++v ) {
            more_ninfo[ v ].dist = f_plus_infinity;
            more_ninfo[ v ].not_included = 1;
        }

        nodeq.init();
        nodeq.put( SOURCE );
        nodeq.put( END_PHASE );

        bool found = true;
        int nphase = 0;

        while ( nphase < n ) {
            int u = nodeq.get();

            if ( END_PHASE == u ) {
                nphase++;

                if ( nodeq.is_empty() ) {
                    found = false;
                    break;
                }

                nodeq.put( END_PHASE );
                continue;
            }
            else 
                more_ninfo[u].not_included = 1;

            float reduced_dist = more_ninfo[ u ].dist - lambda;

            for ( int i = 0; i < g->outdegree( u ); ++i ) {
                int v = g->ith_target_node( u, i );

                float new_dist = reduced_dist + g->ith_target_edge_info( u, i );
                if ( new_dist < more_ninfo[ v ].dist ) {
                    more_ninfo[ v ].dist = new_dist;
                    if ( more_ninfo[ v ].not_included ) {
                        more_ninfo[ v ].not_included = 0;
                        nodeq.put( v );
                    }
                }
            } // for

        } // while nphase > n

        if ( found ) {
            if ( ( upper - lambda ) < EPSILON2 )
                break;
            upper = lambda;
        }
        else {
            if ( ( lambda - lower ) < EPSILON2 )
                break;
            lower = lambda;
        }

    }  // while 

#undef END_PHASE

    delete [] more_ninfo;

    return lambda;
}  // find_min_cycle_mean_for_scc

// End of file
