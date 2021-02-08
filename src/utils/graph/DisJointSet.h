#ifndef __ALGORITHM_DISJOIN_SET_DISJOIN_SET_H__
#define __ALGORITHM_DISJOIN_SET_DISJOIN_SET_H__

#include <vector>
#include <cassert>
#include <boost/config.hpp>
#include <iostream>
#include <algorithm>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
/**********************************************************
 *
 * @Brief :
 *      Get connection component by BGL.
 *      I keep the old function names to suit the old logic codes.
 * *******************************************************/
using namespace boost;
namespace BGIQD{
namespace GRAPH
{

        class DisJoin_Set{
            public :
                void AddConnect( const int a , const int b )
                {
                    add_edge(a, b, G);
                }

                int GetGroup( const int a ) const
                {
                    return results.at(a);
                }

                // call this after add all connections
                // and before any GetGroups
                void GenAllResult() {
                    results.resize(num_vertices(G));
                    connected_components(G, &results[0]);
                }
            private:
                std::vector<int> results;
                typedef adjacency_list< vecS, vecS, undirectedS > Graph;
                Graph G;
        };
} // namespace Algorithm
} // namespace BGIQD
#endif //__ALGORITHM_DISJOIN_SET_DISJOIN_SET_H__
