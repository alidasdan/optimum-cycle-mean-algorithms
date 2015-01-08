//
// COPYRIGHT C 1999- Ali Dasdan (ali_dasdan@yahoo.com)
//
//
// An implementation of Karp's minimum mean cycle algorithm.
//
// @article{Ka78,
// author = {R.~M. Karp},
// title = {A Characterization of the Minimum Cycle Mean in a Digraph},
// journal = {Discrete Mathematics},
// year = 1978,
// volume = 23,
// number = "",
// pages = {309--311},
// month = "",
// }
//

#include "ad_graph.h"

/* ARGSUSED2 */
float
find_min_cycle_mean_for_scc( const ad_graph< ninfo > *g, 
                             int plus_infinity,
                             float lambda_so_far )
{
    int n = g->num_nodes();

    int   *dist_table = new int[ n * ( n + 1 ) ];
    float *max_table = new float[ n ];

    // STEP: Initialize the distances.

    // Set to 0 the distance of length 0 to the source node 0.
    *dist_table = 0;   // D[ 0 ][ 0 ] = 0

    // Initialize D[ k ][ v ] = +infinity for all k and v.
    init_table( dist_table, 1, n * ( n + 1 ) - 1, plus_infinity );

    // STEP: Find the distance of length k from the source node 0 to
    // each node v, for all k.
    int *Dkv = dist_table + n;
    for ( int k = 1; k <= n; ++k ) {
        int *Dk_1 = Dkv - n;  // Dk_1[ v ] = D[ k - 1 ][ v ]

        for ( int v = 0; v < n; ++v ) {
            int *ptr = Dkv++;

            for ( int i = 0; i < g->indegree( v ); ++i )
                min2( *ptr, Dk_1[ g->ith_source_node( v, i ) ] + g->ith_source_edge_info( v, i ) );
        }
    }

    // STEP: Compute lambda using Karp's theorem.
    // lambda = min_{v} max_{0<=k<=n-1} (D_n(k) - Dk(v)) / (n - k).

    // Initialize max_table to -infinity.
    init_table( max_table, 0, n - 1, ( float ) -plus_infinity );

    // Compute the entries of max_table.
    Dkv = dist_table;
    int *Dn = Dkv + n * n;
    int n_minus_k = n;
    for ( int k = 0; k < n; ++k, --n_minus_k )
        for ( int v = 0; v < n; ++v )
            max2( max_table[ v ], ( ( float ) ( Dn[ v ] - *Dkv++ ) ) / n_minus_k );

    // Compute lambda using max_table.
    float lambda = compute_min( max_table, 0, n - 1, ( float ) plus_infinity );

    delete [] dist_table;
    delete [] max_table;

    return lambda;
}  // find_min_cycle_mean_for_scc

// End of file
