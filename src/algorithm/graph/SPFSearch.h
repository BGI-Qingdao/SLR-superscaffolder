#ifndef __ALGORITHM_GRAPH_SPFSEARCH_H__
#define __ALGORITHM_GRAPH_SPFSEARCH_H__
#include <map>
#include "algorithm/graph/Graph.h"
#include "algorithm/fibheap/fib_heap.h"
namespace BGIQD{
    namespace GRAPH{

        template<class BaseNode>
            struct SPFNode  : public BGIQD::FIBHEAP::Node<int,typename BaseNode::NodeNodeId>
            {
                typedef typename BaseNode::NodeNodeId NodeId;
                typedef BGIQD::FIBHEAP::Node<int,typename BaseNode::NodeNodeId> Base;
                enum Type
                {
                    Unknow = 0 ,
                    InHeap = 1 ,
                    Finish = 2 ,
                };
                Type type;

                void InitAsRoot(const BaseNode & node )
                {
                    Base::value = node.id ;
                    Base::key = 0;
                    type = Type::Finish ;
                }
                void Init(const BaseNode & node ,const SPFNode & father , int length)
                {
                    Base::value = node.id ;
                    Base::key = father.key + length;
                    type = Type::InHeap ;
                }

                NodeId prev ;
            };

        template<class GraphAccess
            , class EdgeItr
            , class PathEnder
            , class SPFNode1 = SPFNode<typename GraphAccess::Node > >
            struct SPFSearch 
            {
                typedef typename GraphAccess::GraphNodeId   NodeId;
                typedef typename GraphAccess::GraphEdgeId   EdgeId;
                typedef typename GraphAccess::Node          NodeBase;
                typedef typename GraphAccess::Edge          EdgeBase;

                typedef SPFNode1                            FibNode;
                typedef BGIQD::FIBHEAP::FibHeap<FibNode>    FibHeap;

                void DoSPFSearch(NodeId start)
                {
                    NodeBase & root = accesser.AccessNode(start);

                    //path.push(EdgeItr(accesser.AccessEdge(root.edge_id , root.id) , accesser));
                    //ender.AddEdge(accesser.AccessEdge(root.edge_id));
                    ender.Start();
                    NodeId prev ;
                    fib_nodes[start].InitAsRoot(root);
                    heap.Insert(fib_nodes[start]);
                    while ( ! heap.Empty() )
                    {
                        auto & min = heap.ExtractMin() ;
                        EdgeItr edge_itr(accesser.AccessEdge(root.edge_id , root.id) , accesser);
                        ender.AddEdge(edge_itr);
                        if( ender.IsEnd() )
                        {
                            ender.PopEdge();
                            continue ;
                        }
                        AddEdge(edge_itr);
                        ender.PopEdge();
                    }
                }

                PathEnder ender;

                GraphAccess accesser;

                protected:
                    FibHeap heap;
                    std::map<NodeId , FibNode> fib_nodes;

                    void AddEdge(EdgeItr & itr)
                    {
                        NodeId from = itr->from;
                        FibNode * f_node = NULL ;
                        try{
                            f_node = &fib_nodes.at(from);
                        }
                        catch ( ... )
                        {
                            assert(0);
                        }
                        assert(f_node);

                        f_node->type = FibNode::Type::Finish ;

                        while( itr != EdgeItr::end() )
                        {
                            NodeId to = itr->to ;
                            auto & base_to_node = accesser.AccessNode(to);
                            int l = -1;
                            accesser.GetAttr(base_to_node,*itr,"length",l);
                            assert( l >= 0 );

                            auto itr_n = fib_nodes.find( to );
                            if( itr_n == fib_nodes.end() )
                            {
                                // This is a new node ;
                                auto & new_node = fib_nodes[to];
                                new_node.Init(base_to_node,*f_node);
                                ender.AddNode(base_to_node ,new_node);
                                if( ender.IsEnd() )
                                {
                                    ender.PopNode() ;
                                    continue ;
                                }
                                heap.Insert(new_node);
                                // pop it because push history is useless .
                                ender.PopNode();
                            }
                            else
                            {
                                assert( itr_n -> type == FibNode::Type::InHeap );
                                if( f_node->key + l < itr_n->key )
                                {
                                    heap.DecreaseKey(*itr , f_node->key+l);
                                }
                            }
                        }
                    }
            };
    }
}
#endif //__ALGORITHM_GRAPH_SPFSEARCH_H__
