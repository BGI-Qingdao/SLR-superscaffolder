#ifndef __SOAP2_CONTIGGRAPHSPF_H__
#define __SOAP2_CONTIGGRAPHSPF_H__

#include "soap2/contigGraph.h"
#include "soap2/contigGraphSearch.h"
#include "algorithm/graph/Graph.h"
#include "algorithm/graph/SPFSearch.h"

namespace BGIQD{
    namespace SOAP2{

        struct SFPEnder:
            public BGIQD::GRAPH::PathEndHelperBase<
            GraphEA_Access 
            , traits_search_node
            , BGIQD::GRAPH::SPFNode<typename GraphEA_Access::Node> 
            >
        {
            bool ender_flag ;

            int max_branch ;

            typedef std::function<NodeType(NodeId)> NodeTypeDetector;

            NodeTypeDetector keyer ;

            int max_length ;

            std::map<NodeId , Tos> founder ;

            void AddNode(const Node & node , const SNode & dnode ) 
            {
                if( int(founder.size()) >= max_branch )
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
                    if ( max_length != -1 && dnode.Base.key > max_length )
                    {
                        ender_flag = true ;
                        return ;
                    }
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
                    assert(0);
                ender_flag = true ;
            }

            void PopEdge() { ender_flag  = false ; }

            void PopNode() { ender_flag  = false ; }

            void AddEdge(const Edge & ) { ender_flag = false ; }

            void Start()
            {
                ender_flag = false ;
            }
            void Init( NodeTypeDetector  k , int max_l ,int max_b )
            {
                keyer= k ;
                max_length = max_l;
                max_branch = max_b ;
            }
            bool IsEnd() const { return ender_flag ; }
        };

    }
}

#endif //__SOAP2_CONTIGGRAPHSPF_H__
