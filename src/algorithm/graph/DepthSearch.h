#ifndef __ALGORITHM_GRAPH_DEPTHSEARCH_H__
#define __ALGORITHM_GRAPH_DEPTHSEARCH_H__ 

#include <map>
#include <set>
#include <stack>
#include <sstream>
#include <string>
#include <iostream>
#include <cassert>
namespace BGIQD{
    namespace GRAPH{

        //        template<class Node>
        //            class NodeToIterator
        //            {
        //                public:
        //                    void 
        //
        //            };

        enum DepthSearchEdgeType
        {
            Invalid = -1 ,
            White = 0 ,
            Gray = 1 ,
            Black = 2 ,
            EndPoint  = 3 ,
        };

        template<class NodeBase>
            struct DepthSearchNode
            {
                typedef NodeBase                      Node;
                typedef typename NodeBase::NodeNodeId NodeId;
                typedef DepthSearchEdgeType           Type;
                Type                                type ;

                NodeId                              id;
                NodeId                              prev ;

                int                                 first_found ;
                int                                 finish_search ;

                static const  NodeId                invalid = -1;

                std::set<NodeId>                    backword_from ;
                std::set<NodeId>                    forword_from ;
                std::set<NodeId>                    cross_from ;

                void Clean()
                {
                    type = White ;
                    prev = invalid ;
                    first_found = -1 ;
                    finish_search = -1 ;

                    backword_from.clear();
                    forword_from.clear();
                    cross_from.clear();

                }

                void ReSetParent(const Node & me,const DepthSearchNode & parenet,int step_start )
                {
                    prev = parenet.id ;
                    first_found = step_start ;
                    type = White ;
                    id = me.id ;
                    //TODO : reset edges
                }

                void Init(const Node & me,const DepthSearchNode & parenet,int step_start )
                {
                    type = White ;
                    prev = parenet.id;
                    id = me.id;
                    first_found =  step_start ;
                }

                std::string ToString() const {
                    std::ostringstream ost;
                    ost<<"id = "<<id<<" prev= "<<prev;
                    ost<<" first= "<<first_found<<" last= "<<finish_search;
                    for( auto x : backword_from)
                    {
                        ost<<" Backword= "<<x;
                    }
                    for( auto x : forword_from)
                    {
                        ost<<" Forword= "<<x;
                    }
                    for( auto x : cross_from)
                    {
                        ost<<" Cross= "<<x;
                    }
                    return ost.str();
                }
            };

        template<class GraphAccess , class traits
            , class DepthNode = DepthSearchNode<typename GraphAccess::Node> >
            struct DepthSearchPathEndHelperBase
            {
                typedef typename GraphAccess::GraphNodeId NodeId;
                typedef typename GraphAccess::GraphEdgeId EdgeId;
                typedef typename GraphAccess::Node        Node;
                typedef typename GraphAccess::Edge        Edge;
                typedef DepthNode                         DNode;
                typedef GraphAccess                       Access;
                typedef traits                            traisId;

                void Start() {assert(0);}
                void AddNode(const Node & , const DepthNode &) {assert(0);}
                void AddEdge(const Edge & ) {assert(0);}
                void PopEdge() {assert(0);}
                void PopNode() {assert(0);}
                bool IsEnd() const { assert(0) ;} ;
            };

        template<class GraphAccess
            , class EdgeItr
            , class PathEnder
            , class DepthNode = DepthSearchNode<typename GraphAccess::Node > >
            struct DepthSearch
            {
                typedef typename GraphAccess::GraphNodeId   NodeId;
                typedef typename GraphAccess::GraphEdgeId   EdgeId;
                typedef typename GraphAccess::Node          NodeBase;
                typedef typename GraphAccess::Edge          EdgeBase;
                typedef DepthNode                           Node;
                //              typedef DepthSearchEdge<EdgeBase>           Edge;

                NodeId                                      start;

                std::map<NodeId , Node>                     nodes;
                //std::map<EdgeId , Edge>                   edges;
                std::stack<EdgeItr>                         path;
                GraphAccess                                 accesser;
                PathEnder                                   ender;

                void PrintNodes() const
                {
                    for( auto & i : nodes )
                    {
                        std::cerr<<i.second.ToString()<<std::endl;
                    }
                }
                //  do depth search.
                //      1. stop go deeper if ender say yes.
                //      2. avoid to use recursion because the depth may very large.
                //
                // @params
                //          start : id of root node
                //          step  : index of start step
                // @return 
                //          end step index
                //
                int DoDepthSearch(NodeId start , int s_step )
                {
                    int step = s_step;
                    NodeBase & root = accesser.AccessNode(start);
                    Node & curr = nodes[start];
                    //curr.base = root ;
                    curr.id = start ;
                    curr.first_found = step ; // 0 as root ;
                    curr.type = Node::Type::White;
                    curr.prev = Node::invalid ;

                    path.push(EdgeItr(accesser.AccessEdge(root.edge_id , root.id) , accesser));
                    //ender.AddEdge(accesser.AccessEdge(root.edge_id));
                    ender.Start();
                    bool new_node_in_path = true ;
                    NodeId prev ;
                    while ( ! path.empty() )
                    {
                        auto & itr = path.top() ;
                        Node & top = nodes.at( itr.node_id) ;

                        if( itr == EdgeItr::end() )
                        {
                            step ++ ;
                            if( top.type == Node::Type::White )
                                top.type = Node::Type::EndPoint ;
                            else
                                top.type = Node::Type::Black ;
                            top.finish_search = step ;
                            path.pop() ;
                            new_node_in_path = false ;
                            ender.PopEdge();
                            ender.PopNode();
                            continue ;
                        }

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

                            ender.AddNode(accesser.AccessNode(top.id),top);
                            if( ender.IsEnd() )
                            {
                                if( top.type == Node::Type::White )
                                { // end by new node
                                    top.type = Node::Type::EndPoint ;
                                    step ++ ;
                                    top.finish_search = step ;
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

                            if( top.type != Node::Type::White )
                            {
                                step ++ ;
                                top.ReSetParent( accesser.AccessNode(top.id) , nodes.at(prev), step );
                            }
                            //in path now , force make it Gray !!!
                            top.type = Node::Type::Gray ;
                        }

                        NodeId next_node_id  = itr->to;
                        auto itr_n = nodes.find( next_node_id ) ;
                        if ( itr_n == nodes.end() || itr_n->second.type == Node::Type::White )
                        {
                            step ++ ;
                            // white node alloc 
                            auto & next_node = nodes[next_node_id] ;
                            //next_node.id = next_node_id;
                            //next_node.type = Node::Type::White ;
                            //next_node.first_found = step ;
                            //next_node.prev = top.id ;

                            next_node.Init( accesser.AccessNode(next_node_id) , top, step );
                        }
                        else
                        {
                            auto & node = itr_n->second ;
                            if ( node.type == Node::Type::Gray )
                            { // match a backwoard edge
                                node.backword_from.insert( top.id ) ;
                            }
                            else
                            {// black node or endpoint left

                                if( node.first_found > top.first_found  )
                                { // to child ,match a forward edge
                                    node.forword_from.insert(top.id);
                                }
                                else
                                { // to other branch ,match a corss edge
                                    node.cross_from.insert(top.id);
                                }
                            }
                        }
                        prev = top.id ;
                        auto & next_node = accesser.AccessNode( next_node_id) ;
                        EdgeId next_edge_id = next_node.edge_id;
                        EdgeBase & next_edge = accesser.AccessEdge( next_edge_id ,next_node_id) ;

                        path.push(EdgeItr(next_edge, accesser));
                        ender.AddEdge(next_edge);

                        new_node_in_path = true ;
                        ++ itr ;
                    }
                    return step ;
                }
            };
    } // namespace GRAPH
} // namespace BGIQD

#endif //__ALGORITHM_GRAPH_DEPTHSEARCH_H__
