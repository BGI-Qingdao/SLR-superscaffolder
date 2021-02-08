#ifndef __ALGORITHM_GRAPH_MINTREE_H__
#define __ALGORITHM_GRAPH_MINTREE_H__

#include <boost/config.hpp>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <cassert>
/**********************************************************
 *
 * @Brief :
 *
 *  call prim's mininum spanning tree algorithm in boost 
 *  graph library.
 *
 * *******************************************************/
using namespace boost;
namespace BGIQD {
    namespace GRAPH {

        // User implement this to get weight:
        template<class TNode,class TEdge, class Value>
            struct NodeAttr
            {
                Value GetValue( const TEdge & ) const {}
            };


        template<class TListGraph 
            , class TValue 
            , class TNodeAttr =
            NodeAttr<typename TListGraph::Node 
            ,typename TListGraph::Edge
            , TValue> >
            struct MinTreeHelper
            {
                typedef TListGraph  TheGraph;
                typedef TNodeAttr   TheAttr;
                typedef TValue      Value;

                typedef typename TListGraph::Node Node;
                typedef typename TListGraph::Edge Edge;
                typedef typename Edge::EdgeEdgeId EdgeId;
                typedef typename Edge::EdgeNodeId NodeId;

                // Call MST in BGL
                TheGraph MinTree(const TheGraph & base , const TheAttr & attr )
                {
                    TheGraph ret;
                    typedef adjacency_list< vecS, vecS, undirectedS,
                            property< vertex_distance_t, int >, property< edge_weight_t, float> >
                                BGLGraph;
                    BGLGraph G;
                    // Init a boost graph
                    for(const auto & edge : base.edges ){
                        if(  edge.IsValid() ){
                            auto u = vertex(edge.from,G);
                            auto v = vertex(edge.to,G);
                            add_edge(u,v,attr.GetValue(edge),G);
                        }
                    }
                    // call mst
                    std::vector< graph_traits<BGLGraph >::vertex_descriptor > p(num_vertices(G));
                    prim_minimum_spanning_tree(G, &p[0]);

                    // construct the return graph to fit the old codes
                    std::set<EdgeId> edge_to_del;
                    for(const auto & edge : base.edges ){
                        if(p[edge.from] == edge.to || p[edge.to] == edge.from )
                            continue;
                        edge_to_del.insert(edge.id);
                    }
                    ret=base;
                    for( EdgeId eid: edge_to_del)
                    {
                        ret.RemoveEdge(eid);
                    }
                    return ret ;
                }
            };
    }
}

#endif //__ALGORITHM_GRAPH_MINTREE_H__
