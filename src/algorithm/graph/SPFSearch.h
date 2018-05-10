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
                typename BaseNode::NodeNodeId NodeId;
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
                    bool new_node_in_path = true ;
                    NodeId prev ;
                    while ( ! heap.Empty() )
                    {

                    }
                }
                PathEnder ender;
                protected:
                    FibHeap heap;
                    std::map<NodeId , FibNode> fib_nodes;

                    void AddEdge(EdgeItr & itr)
                    {
                        while( itr != EdgeItr::end() )
                        {
                            NodeId to = itr->to ;
                            auto itr = fib_nodes.find( to );
                            if( itr == fib_nodes.end() )
                            {

                            }
                        }
                    }
            };
    }
}
#endif //__ALGORITHM_GRAPH_SPFSEARCH_H__
