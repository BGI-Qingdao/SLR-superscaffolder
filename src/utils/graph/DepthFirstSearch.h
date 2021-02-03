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

        template<class Graph>
            struct EdgeIterator
        {
            public:

                typedef typename Graph::EdgeId        Id;
                typedef typename Graph::Edge          Edge;

                EdgeIterator() : curr(NULL) ,accessor( NULL ) { }

                EdgeIterator(const Edge & e , Graph & acc )
                {
                    node_id = e.from;
                    if( e.id != Edge::invalid )
                        curr = &e;
                    else
                        curr = NULL ;
                    accessor = &acc ;
                }

                EdgeIterator( const EdgeIterator & ei )
                {
                    node_id = ei.node_id;
                    curr = ei.curr ;
                    accessor = ei.accessor ;
                }

                EdgeIterator & operator = ( const EdgeIterator & ei )
                {
                    if( &ei != this )
                    {
                        node_id = ei.node_id ;
                        curr = ei.curr ;
                        accessor = ei.accessor ;
                    }
                    return *this;
                }

                // ONLY detect curr.
                // ONLY use == with End() .
                bool operator == ( const EdgeIterator & ei )const
                {
                    return curr == ei.curr ;
                }

                // ONLY detect curr.
                // ONLY use == with End() .
                bool operator != ( const EdgeIterator & ei )const
                {
                    return curr != ei.curr ;
                }

                const Edge & operator*() const  { return *curr ; }

                const Edge * operator->() const  { return curr ; }

                EdgeIterator & operator ++() {
                    if( curr != NULL && accessor != NULL )
                    {
                        Id next = curr->next ;
                        if( next != Edge::invalid )
                            curr = &(accessor->GetEdge(next));
                        else
                            curr = NULL ;
                    }
                    else
                        curr = NULL ;
                    return *this ;
                }

                static EdgeIterator & end() 
                {
                    static EdgeIterator end;
                    return end ;
                }
                typename Edge::EdgeNodeId node_id ;
            private:
                const Edge * curr ;
                Graph      * accessor ;
        };

        enum DepthSearchEdgeType
        {
            Invalid = -1 ,
            White = 0 ,
            Gray = 1 ,
            Black = 2 ,
            EndPoint  = 3 ,
        };

        template<class Node>
            struct DepthSearchNode
            {
                typedef typename Node::NodeNodeId   NodeId;
                typedef DepthSearchEdgeType         Type;
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

                void InitRoot(const Node & me,int step_start )
                {
                    type = White ;
                    prev = invalid;
                    id = me.id;
                    first_found =  step_start ;
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

        template<class Graph
            , class EdgeItr
            , class PathEnder
            , class DepthNode = DepthSearchNode<typename Graph::Node > >
            struct DepthSearch
            {
                typedef typename Graph::NodeId   NodeId;
                typedef typename Graph::EdgeId   EdgeId;
                typedef typename Graph::Node          NodeBase;
                typedef typename Graph::Edge          EdgeBase;
                typedef DepthNode                           Node;
                //              typedef DepthSearchEdge<EdgeBase>           Edge;

                NodeId                                      start;

                std::map<NodeId , Node>                     nodes;
                //std::map<EdgeId , Edge>                   edges;
                std::stack<EdgeItr>                         path;
                Graph                                       *accesser;
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
                    NodeBase & root = accesser->GetNode(start);
                    Node & curr = nodes[start];

                    curr.InitRoot(root,step);

                    path.push(EdgeItr(accesser->GetEdge(root.edge_id ) , *accesser));
                    ender.Start();
                    bool new_node_in_path = true ;
                    NodeId prev ;
                    while ( ! path.empty() )
                    {
                        auto & itr = path.top() ;
                        Node * test ;
                        try {
                            Node & top = nodes.at( itr.node_id ) ;
                            test = &top;
                        }
                        catch(...)
                        {
                            assert(0);
                        }
                        Node & top = *test;
                        if( itr == EdgeItr::end() )
                        {
                            step ++ ;
                            if( top.type == Node::Type::White )
                                top.type = Node::Type::EndPoint ;
                            else
                                top.type = Node::Type::Black ;
                            top.finish_search = step ;
                            path.pop() ;
                            ender.PopEdge();
                            if( ! new_node_in_path )
                                ender.PopNode();

                            new_node_in_path = false ;
                            //assert(path.size() == ender.nodes.size() );
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
                                //assert(path.size() == ender.nodes.size() );
                                continue ;
                            }

                            ender.AddNode(accesser->GetNode(top.id),top);
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
                                //assert(path.size() == ender.nodes.size() );
                                continue ;
                            }

                            if( top.type != Node::Type::White )
                            {
                                step ++ ;
                                try {
                                    top.ReSetParent( accesser->GetNode(top.id) , nodes.at(prev), step );
                                }
                                catch(...)
                                {
                                    assert(0);
                                }
                            }
                            //in path now , force make it Gray !!!
                            top.type = Node::Type::Gray ;
                        }

                        //assert(path.size() == ender.nodes.size() );
                        NodeId next_node_id  = itr->to;
                        auto itr_n = nodes.find( next_node_id ) ;
                        if ( itr_n == nodes.end() || itr_n->second.type == Node::Type::White )
                        {
                            step ++ ;
                            // white node alloc 
                            auto & next_node = nodes[next_node_id] ;
                            next_node.Init( accesser->GetNode(next_node_id) , top, step );
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
                        auto & next_node = accesser->GetNode( next_node_id) ;
                        EdgeId next_edge_id = next_node.edge_id;
                        EdgeBase & next_edge = accesser->GetEdge( next_edge_id) ;
                        //assert(path.size() == ender.nodes.size() );
                        path.push(EdgeItr(next_edge, *accesser));
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
