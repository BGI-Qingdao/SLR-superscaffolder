#ifndef __SOAP2_CONTIGGRAPHDEPTH_H__
#define __SOAP2_CONTIGGRAPHDEPTH_H__

#include "algorithm/graph/DepthSearch.h"
#include "algorithm/graph/Graph.h"
#include "soap2/contigGraph.h"
namespace BGIQD{
    namespace SOAP2{

        struct Node_EA : public BGIQD::GRAPH::GraphNodeBase<unsigned int , long>
        {
            int length ;
        };

        struct DNode_EA : public BGIQD::GRAPH::DepthSearchNode<Node_EA>
        {
            int path_length ;
        };

        struct GraphEA_Access : public BGIQD::GRAPH::GraphAccessBase<
                                BGIQD::SOAP2::GraphEA
                                , unsigned int  
                                , long
                                , Node_EA
                                >
        {
            Node & AccessNode(GraphNodeId i)
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

            Edge & AccessEdge(GraphEdgeId i , GraphNodeId from)
            {
                static Edge none;
                if( i == Edge::invalid )
                {
                    assert(0);
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
        };

        struct traits_search_node {} ;

        struct traits_search_path {} ;

        struct DepthSearchEAEnder :
            public BGIQD::GRAPH::DepthSearchPathEndHelperBase<
                            GraphEA_Access 
                            , traits_search_node
                            >
        {
            private:
                int curr_length ;
                int curr_depth ;
                bool ender_flag ;

                int max_length ;
                int max_depth ;

                std::stack<Node> nodes;

            public:

                enum NodeType
                {
                    Unknow = 0 ,
                    Normal = 1 ,
                    Key_Unknow = 2,
                    RC_Key_Unknow = 3,
                    Key_Neibs  = 4 , 
                    RC_Key_Neibs = 5 ,
                };

                typedef std::function<NodeType(NodeId)> NodeTypeDetector;

                struct Tos
                {
                    bool base ;
                    bool bal ;
                    Tos() : base(false) , bal( false) {}
                } ;

                std::map<NodeId,Tos> founder;

                void Init( NodeTypeDetector key ,int max_l , int max_d )
                {
                    max_length = max_l ;
                    max_depth = max_d ;
                    keyer = key ;
                }


                void Start() 
                {
                    curr_length = 0 ;
                    curr_depth = 0;
                    ender_flag = false;
                }

                void AddNode(const Node & node , BGIQD::GRAPH::DepthSearchEdgeType ) 
                {
                    if( nodes.empty() )
                    {
                        nodes.push(node);
                        return ;
                    }
                    nodes.push(node) ;
                    curr_depth ++ ;
                    curr_length += node.length ;
                    if ( max_depth != -1 && curr_depth > max_depth )
                    {
                        ender_flag = true ;
                        return ;
                    }
                    if ( max_length != -1 && curr_length > max_length )
                    {
                        ender_flag = true ;
                        return ;
                    }

                    auto ret =  keyer(node.id);
                    if( ret == NodeType::Unknow )
                    {
                        assert(0);
                    }
                    if( ret == NodeType::Normal ) 
                    {
                        return ;
                    }
                    else if ( ret == NodeType::Key_Neibs)
                    {
                        founder[node.id].base = true ;
                    }
                    else if ( ret == NodeType::RC_Key_Neibs )
                    {
                        founder[node.id-1].bal = true ;
                    }
                    else if( ret == NodeType::Key_Unknow || ret == NodeType::RC_Key_Unknow )
                    {
                        // Do nothing 
                        ;
                    }
                    else
                        ;
                    ender_flag = true ;
                }

                void AddEdge(const Edge & ) 
                {

                }

                void PopEdge() {

                }
                void PopNode() {
                    ender_flag  =  false ;
                    if( nodes.empty() )
                    {
                        //assert(0);
                    }
                    else
                    {
                        auto & t = nodes.top() ;
                        curr_depth -- ;
                        curr_length -= t.length ;
                        nodes.pop();
                    }
                    assert( curr_depth >= 0 );
                    assert( curr_length>= 0 );
                }

                bool IsEnd() const {
                    return ender_flag ;
                }

            private:
                NodeTypeDetector keyer ;
        };

    } // namespace SOAP2
} // namespace BGIQD

#endif //__SOAP2_CONTIGGRAPHDEPTH_H__
