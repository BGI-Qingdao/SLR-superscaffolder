#ifndef __STLFR_CONTIGPEGRAPHSEARCH_H__
#define __STLFR_CONTIGPEGRAPHSEARCH_H__

#include "utils/graph/DepthFirstSearch.h"
#include "stLFR/contigPEGraph.h"

#include <functional>

namespace BGIQD{
    namespace stLFR{

        struct  ContigPEGraphAccess  : public BGIQD::GRAPH::GraphAccessBase<
                                BGIQD::stLFR::ContigPEGraph
                                , unsigned int
                                , long
                                , BGIQD::stLFR::ContigNode
                                , BGIQD::stLFR::PEEdge
                                >
            {

                Node & AccessNode(const GraphNodeId & id) 
                {
                    return (*base).GetNode(id);
                }

                Edge invalid ;

                void Init()
                {
                    invalid.id = Edge::invalid ;
                }

                Edge & AccessEdge(const GraphEdgeId & id, const GraphNodeId & nodeId)
                {
                    if( id == Edge::invalid )
                    {
                        invalid.from = nodeId ;
                        return invalid ;
                    }
                    return (*base).GetEdge(id);
                }

            };

        struct traits_search_node{};

        struct DepthEnder:
            public BGIQD::GRAPH::PathEndHelperBase<
            ContigPEGraphAccess 
            , traits_search_node
            , BGIQD::GRAPH::DepthSearchNode<typename ContigPEGraphAccess::Node> 
            >
        {

            typedef BGIQD::GRAPH::PathEndHelperBase<
            ContigPEGraphAccess 
            , traits_search_node
            , BGIQD::GRAPH::DepthSearchNode<typename ContigPEGraphAccess::Node> 
            > Basic;
            bool ender_flag ;
            typedef typename Basic::NodeId NodeId;

            enum NodeType
            {
                Unknow = 0 ,
                Others = 1 ,
                Seed = 2 ,
            };

            typedef std::function<NodeType(const NodeId &)> NodeTypeDetector;

            NodeTypeDetector keyer ;

            int max_length ;

            int curr_len ;

            int last_len ;

            bool first_in  ;

            std::stack<const Node * > path;

            void AddNode(const Node & node , const SNode & snode )
            {
                path.push(&node);

                if( ! first_in )
                {
                    first_in = true ;
                    return ;
                }
                if( snode.type != BGIQD::GRAPH::DepthSearchEdgeType::White )
                {
                    ender_flag = true ;
                    return ;
                }
                auto ret =  keyer(node.id);
                if( ret == NodeType::Unknow )
                {
                    assert(0);
                }
                if( ret == NodeType::Others) 
                {
                    if ( max_length != -1 && curr_len > max_length )
                    {
                        ender_flag = true ;
                        return ;
                    }
                    else
                    {
                        curr_len += node.contigLen;
                        return ;
                    }
                }
                else if ( ret == NodeType::Seed)
                {
                    ender_flag = true ;
                    return ;
                }
                else
                    assert(0);
                ender_flag = true ;
            }

            void PopEdge() { ender_flag  = false ; }

            void PopNode()
            {
                if( path.empty() )
                {
                    assert(0);
                    return ;
                }
                const Node * top = path.top();
                path.pop();
                ender_flag  = false ;
                curr_len = top->contigLen ;
            }

            void AddEdge(const Edge & ) { ender_flag = false ; }

            void Start()
            {
                ender_flag = false ;
                curr_len =last_len = 0 ;
            }
            void Init( NodeTypeDetector  k , int max_l )
            {
                keyer= k ;
                max_length = max_l;
                first_in = false ;
            }
            bool IsEnd() const { return ender_flag ; }
        };

    } // namespace SOAP2
} // namespace BGIQD
#endif
