#ifndef __ALGORITHM_GRAPH_DEPTHSEARCH_H__
#define __ALGORITHM_GRAPH_DEPTHSEARCH_H__ 

#include <map>
#include <set>
#include <stack>
namespace BGIQD{
    namespace GRAPH{

        //        template<class Node>
        //            class NodeToIterator
        //            {
        //                public:
        //                    void 
        //
        //            };

        template<class NodeBase>
            struct DepthSearchNode
            {
                typedef typename NodeBase::NodeNodeId NodeId;
                enum Type
                {
                    White = 0 ,
                    Gray = 1 ,
                    Black = 2 ,
                    EndPoint  = 3 ,
                };
                Type                                type ;

                NodeId                              id;
                NodeId                              prev ;

                int                                 first_found ;
                int                                 finish_search ;

                static NodeId                       invalid = -1;

                std::set<NodeId>                    backword_from ;
                std::set<NodeId>                    forword_from ;
                std::set<NodeId>                    cross_from ;

                void Clean()
                {
                    type = White ;
                    prev = invalid ;
                    first_found = -1 ;
                    finish_search = -1 ;
                }
            };

        template<class EdgeBase>
            struct DepthSearchEdge
            {
                enum Type
                {
                    TreeEdge = 0 ,
                    Backward = 1 ,
                    Forward = 2 ,
                    Crossward = 3 ,
                };
                Type        type ;

                EdgeBase    base ;
            };

        template<class GraphAccess>
            struct DepthSearchPathEndHelper
            {
                typedef typename GraphAccess::GraphNodeId NodeId;
                typedef typename GraphAccess::GraphEdgeId EdgeId;
                typedef typename GraphAccess::Node        Node;
                typedef typename GraphAccess::Edge        Edge;

                void AddNode(const Node & ) {} 
                void AddEdge(const Edge & ) {}

                void PopEdge() {}
                void PopNode() {}

                bool IsEnd() const { return true ; }
            };

        template<class GraphAccess ,class EdgeItr>
            struct DepthSearch
            {
                typedef typename GraphAccess::GraphNodeId   NodeId;
                typedef typename GraphAccess::GraphEdgeId   EdgeId;
                typedef typename GraphAccess::Node          NodeBase;
                typedef typename GraphAccess::Edge          EdgeBase;
                typedef DepthSearchNode<NodeBase>           Node;
                typedef DepthSearchEdge<EdgeBase>           Edge;

                NodeId                                    start;

                std::map<NodeId , Node>                   nodes;
                //std::map<EdgeId , Edge>                   edges;
                std::stack<EdgeItr>                        path;
                GraphAccess                               accesser;
                DepthSearchPathEndHelper<GraphAccess>     ender;

                //  do depth search.
                //      1. stop go deeper if ender say yes.
                //      2. avoid to use recursion because the depth may very large.
                void DoDepthSearch()
                {
                    int step = 0 ;

                    NodeBase root = accesser.AccessNode(start);
                    Node & curr = nodes[start];
                    //curr.base = root ;
                    curr.id = start ;
                    curr.first_found = step ; // 0 as root ;
                    curr.type = Node::Type::White;
                    curr.prev = Node::invalid ;

                    path.push(EdgeItr(accesser.AccessEdge(root.edge_id) , accesser));

                    bool new_node_in_path = true ;
                    while ( ! path.empty() )
                    {
                        step ++ ;
                        auto & itr = path.top() ;
                        Node & top = nodes.at( itr->from ) ;

                        if( new_node_in_path )
                        {
                            if( ender.IsEnd() )
                            { // end by edge
                                top.Clean();
                                ender.PopEdge();
                                path.pop() ;
                                new_node_in_path = false ;
                                continue ;
                            }
                            ender.AddNode(accesser.AccessNode(top.id));
                            if( ender.IsEnd() )
                            {
                                if( top.type == Node::Type::White )
                                { // end by new node
                                    top.type = Node::Type::EndPoint ;
                                }
                                else 
                                { // end by old node
                                    ;
                                }
                                ender.PopEdge();
                                ender.PopNode();
                                path.pop() ;
                                new_node_in_path = false ;
                                continue ;
                            }
                            //in path now , force make it Gray !!!
                            top.type = Node::Type::Gray ;
                        }
                        if( itr == EdgeItr::end )
                        {
                            if( top.type == Node::Type::White )
                                top.type = Node::Type::EndPoint ;
                            else
                                top.type = Node::Type::Black ;
                            top.finish_search = step ;
                            path.pop() ;
                            new_node_in_path = false ;
                            continue ;
                        }

                        NodeId curr_node_id  = top.base.edge.to ;
                        EdgeBase curr_edge = accesser.AccessEdge( top.base.edge.id ) ;
                        auto itr_n = nodes.find( curr_node_id ) ;
                        if ( itr_n == nodes.end() || itr->type == Node::Type::White )
                        {
                            // white node alloc 
                            auto & next_node = nodes[curr_node_id] ;
                            next_node.base = accesser.AccessNode(curr_node_id);
                            next_node.type = Node::Type::White ;
                            next_node.first_found = step ;
                            next_node.prev = top.id ;

                        }
                        else
                        {
                            if ( itr->type == Node::Type::Gray )
                            { // match a backwoard edge
                                itr->backword_from.insert( top.id ) ;
                            }
                            else
                            {// black node or endpoint left

                                if( itr->first_found > top.first_found  )
                                { // to child ,match a forward edge
                                    itr->forword_from.insert(top.id);
                                }
                                else
                                { // to other branch ,match a corss edge
                                    itr->cross_form.insert(top.id);
                                }
                            }
                        }
                        path.push(EdgeItr(accesser.AccessEdge(curr_node_id) , accesser));
                        //path.push(curr_node_id);
                        ender.AddEdge(curr_edge);
                        new_node_in_path = true ;
                        itr ++ ;
                    }
                }
            };
    } // namespace GRAPH
} // namespace BGIQD

#endif //__ALGORITHM_GRAPH_DEPTHSEARCH_H__
