#ifndef __SOAP2_CONTIGGRAPHDEPTH_H__
#define __SOAP2_CONTIGGRAPHDEPTH_H__

#include "algorithm/graph/DepthSearch.h"
#include "algorithm/graph/Graph.h"
#include "soap2/contigGraph.h"
#include "soap2/contigGraphSearch.h"

namespace BGIQD{
    namespace SOAP2{

        struct DNode_EA : public BGIQD::GRAPH::DepthSearchNode<Node_EA>
        {
            int path_length ;
            void ReSetParent(const Node & me,const DNode_EA & parenet,int step_start )
            {
                BGIQD::GRAPH::DepthSearchNode<Node_EA>::ReSetParent(me,parenet ,step_start);
                if( parenet.id == invalid )
                {
                    path_length = 0 ;
                }
                else
                {
                    path_length = parenet.path_length + me.length ;
                }
            }

            void InitRoot(const Node & me,int step_start )
            {
                BGIQD::GRAPH::DepthSearchNode<Node_EA>::InitRoot(me,step_start);
                path_length = 0 ;
            }

            void Init(const Node & me,const DNode_EA & parenet,int step_start )
            {
                BGIQD::GRAPH::DepthSearchNode<Node_EA>::Init(me,parenet ,step_start);
                if( parenet.id == invalid )
                {
                    path_length = 0 ;
                }
                else
                {
                    path_length = parenet.path_length + me.length ;
                }
            }
        };


        struct DepthSearchEAEnder :
            public BGIQD::GRAPH::PathEndHelperBase<
                            GraphEA_Access 
                            , traits_search_node
                            , DNode_EA
                            >
        {
            private:

                int curr_length ;


                bool ender_flag ;

                int max_length ;

                int max_depth ;


            public:

                std::stack<SNode> nodes;

                int curr_depth ;


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

                void AddNode(const Node & node , const SNode & dnode ) 
                {
                    nodes.push(dnode) ;
                    curr_depth ++ ;
                    if( nodes.size() == 1)
                    {
                        curr_length = 0 ;
                        return ;
                    }

                    if( dnode.type != BGIQD::GRAPH::DepthSearchEdgeType::White )
                    {
                        // check if this is a shorter path
                        if( dnode.path_length  > nodes.top().path_length - node.length )
                        {
                            ;
                        }
                        else
                        {
                            ender_flag = true ;
                        }
                    }

                    curr_length = dnode.path_length;

                    auto ret =  keyer(node.id);
                    if( ret == NodeType::Unknow )
                    {
                        assert(0);
                    }
                    if( ret == NodeType::Normal ) 
                    {
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
                        curr_depth -- ;
                        nodes.pop();
                        if(! nodes.empty() )
                            curr_length = nodes.top().path_length ;
                        else
                            curr_length = 0 ;
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
