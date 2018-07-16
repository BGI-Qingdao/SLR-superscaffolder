#ifndef __SOAP2_CONTIGGRAPHSEARCH_H__
#define __SOAP2_CONTIGGRAPHSEARCH_H__

#include "algorithm/graph/Graph.h"
#include "soap2/contigGraph.h"

namespace BGIQD {
    namespace SOAP2 {

        struct Node_EA : public BGIQD::GRAPH::GraphNodeBase<unsigned int , long>
        {
            int length ;
            enum EndType 
            {
                Unknow = 0 ,
                KeyEnd = 1 ,
                MaxLength = 2 ,
                MaxDepth = 3 ,
            } ;
            EndType type ;
        };

        struct traits_search_node {} ;

        struct traits_search_path {} ;

        struct GraphEA_Access : public BGIQD::GRAPH::GraphAccessBase<
                                BGIQD::SOAP2::GraphEA
                                , unsigned int  
                                , long
                                , Node_EA
                                >
        {
            int K ;
            Node & AccessNode(const GraphNodeId &i)
            {
                auto itr = nodes.find(i); 
                if( itr == nodes.end() )
                {
                    auto & n = nodes[i] ;
                    auto & e = (*base).edge_array[i] ;
                    n.id = e.id ;
                    assert( e.id >= 0 );
                    if( e.arc != NULL )
                        n.edge_id= e.arc- base->arc_array;
                    else
                        n.edge_id = Edge::invalid ;
                    n.length = e.length ;
                    return n;
                }
                return itr->second ;
            }

            Edge & AccessEdge(const GraphEdgeId & i , const GraphNodeId & from)
            {
                if( i == Edge::invalid )
                {
                    auto & none = edges_ends[from]; 
                    none.id = Edge::invalid ;
                    none.from = from ; 
                    none.to = 0 ;
                    none.next = Edge::invalid ;
                    return none ;
                }
                auto itr = edges.find(i);
                if( itr == edges.end() )
                {
                    auto & n = edges[i] ;
                    auto & base_edge = base->arc_array[i];
                    n.id = i ;
                    n.from = from;
                    n.to   = base_edge.to;
                    if( base_edge.next != NULL )
                        n.next = base_edge.next - base->arc_array;
                    else
                        n.next = Edge::invalid ;
                    return n;
                }
                return itr->second ;
            }

            void GetAttr(const Node & node, const Edge & ,
                    const std::string & name, int &v ) {
                if(  name == "length" )
                {
                    if( node.length > K )
                        v = node.length - K;
                    else
                        v = 0 ;
                }
                else
                    assert(0);
            }
        };

        enum NodeType
        {
            Unknow = 0 ,
            Normal = 1 ,
            Key_Unknow = 2,
            RC_Key_Unknow = 3,
            Key_Neibs  = 4 , 
            RC_Key_Neibs = 5 ,
        };

        struct Tos
        {
            bool base ;
            bool bal ;
            Tos() : base(false) , bal( false) {}
        } ;

    }
}

#endif //__SOAP2_CONTIGGRAPHSEARCH_H__
