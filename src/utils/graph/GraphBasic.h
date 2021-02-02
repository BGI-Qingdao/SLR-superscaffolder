#ifndef __ALGORITHM_GRAPH_H__
#define __ALGORITHM_GRAPH_H__

#include <vector>
#include <set>
#include <map>
#include <cassert>
#include <list>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>

namespace BGIQD {
    namespace GRAPH {

        template<class NodeId ,  class EdgeId >
            struct IGraphEdgeBasic
            {
                typedef NodeId EdgeNodeId;
                typedef EdgeId EdgeEdgeId;

                EdgeId      id ;

                NodeId      from ;
                NodeId      to ;

                bool operator == (const  IGraphEdgeBasic & i )
                {
                    return ( (from == i.from && to == i.to ) 
                            || ( from == i.to && to == i.from ) ) ;
                }

                bool operator != ( const  IGraphEdgeBasic & i )
                {
                    return ! operator == ( i ) ;
                }

                NodeId OppoNode(const NodeId & one) const 
                {
                    if( one == from ) 
                        return to ;
                    else if ( one == to )
                        return from ;
                    assert(0);
                }

                std::string AttrString() const
                {
                    std::ostringstream ost;
                    ost<<" id = "<<id;
                    return ost.str();
                }
                std::string ToString() const
                {
                    std::ostringstream ost;
                    ost<<from<<"\t--\t"<<to<<" [ "<<AttrString()<<" ]";
                    return ost.str();
                }

                bool InvalidMe()
                {
                    if ( id != invalid ) 
                    {
                        id = invalid ;
                        return true ;
                    }
                    else 
                        return false ;
                }

                bool IsValid() const 
                {
                    return id != invalid ;
                }

                static std::string DOTHead()
                {
                    return "graph {";
                }

                static const EdgeId invalid = -1 ;
            };

        template<class NodeId ,  class EdgeId >
            struct IDigraphEdgeBase : public IGraphEdgeBasic<NodeId, EdgeId>
            {
                typedef IGraphEdgeBasic<NodeId, EdgeId> BaseType ;

                typedef NodeId EdgeNodeId;
                typedef EdgeId EdgeEdgeId;

                bool operator == (const  IDigraphEdgeBase & i )
                {
                    return ( (BaseType::from == i.from && BaseType::to == i.to ));
                }

                bool operator != ( const IDigraphEdgeBase & i )
                {
                    return ! operator == ( i ) ;
                }
                std::string AttrString() const
                {
                    return BaseType::AttrString();
                }
                std::string ToString() const
                {
                    std::ostringstream ost;
                    ost<<BaseType::from<<"\t->\t"<<BaseType::to<<" [ "<<AttrString()<<" ]";
                    return ost.str();
                }
                static std::string DOTHead()
                {
                    return "digraph {";
                }
            };

        template<class NodeId ,  class EdgeId >
            struct IGraphNodeBasic
            {
                typedef NodeId NodeNodeId;
                typedef EdgeId NodeEdgeId;

                NodeId                          id ;
                std::set<EdgeId>               edge_ids;

                void AddEdge(const NodeEdgeId &id )
                {
                    edge_ids.insert(id);
                }

                bool HasEdge(const NodeEdgeId &id )
                {
                    return edge_ids.find(id) != edge_ids.end() ; 
                }

                size_t EdgeNum() const 
                {
                    return edge_ids.size() ;
                }

                void CleanEdges()
                {
                    edge_ids.clear();
                }

                bool RemoveEdge( const NodeEdgeId &id )
                {
                    if ( edge_ids.find(id ) != edge_ids.end () )
                    {
                        edge_ids.erase( id ) ;
                        return true ;
                    }
                    else
                        return false ;
                }
            };

        template<class TNode 
            , class TEdge 
            , class TNodes = std::map<typename TNode::NodeNodeId , TNode>
            , class TEdges = std::vector<TEdge> 
            >
            struct ListGraphBasic
            {
                typedef TNode Node ;
                typedef TEdge Edge ;
                typedef TNodes Nodes ;
                typedef TEdges Edges ;

                typedef typename Node::NodeNodeId NodeId ;
                typedef typename Edge::EdgeEdgeId EdgeId ;

                Nodes nodes ;
                Edges edges ;

                void Init()
                {
                }

                // DO NOT use this if your TNode has more member variable.
                void AddNode(const NodeId &id )
                {
                    nodes[id].id = id ;
                    //nodes[id].edge_id = Edge::invalid ;
                }

                void AddNode(const Node & n )
                {
                    nodes[n.id] = n ;
                }

                Node & GetNode(const NodeId & id )
                {
                    return nodes[id] ;
                }

                const Node & GetNode(const NodeId & id )const 
                {
                    return nodes.at(id) ;
                }

                bool HasNode( const NodeId & id ) const
                {
                    return nodes.find(id) != nodes.end() ;
                }

                Edge & GetEdge( const EdgeId & id )
                {
                    return edges[id] ;
                }

                bool RemoveNode( const NodeId & id )
                {
                    if ( nodes.find(id) != nodes.end() )
                    {
                        auto & n1 = GetNode( id );
                        for( auto  x : n1.edge_ids )
                        {
                            RemoveEdge(x);
                        }
                        nodes.erase(id);
                        return true ;
                    }
                    else
                        return false ;
                }

                bool RemoveEdge( const EdgeId& id )
                {
                    auto & edge = GetEdge( id );
                    if( edge.IsValid() )
                    {
                        auto & n1 = GetNode( edge.from );
                        auto & n2 = GetNode( edge.to );
                        n1.RemoveEdge(id);
                        n2.RemoveEdge(id);
                        edge.InvalidMe();
                        return true ;
                    }
                    else
                        return false ;
                }

                const Edge & GetEdge( const EdgeId &id )const
                {
                    return edges[id] ;
                }
                bool CheckEdge(const NodeId & from ,const  NodeId & to )
                {
                    Edge tmp ;
                    tmp.from = from ;
                    tmp.to = to ;
                    if( ! HasNode( from ) || !HasNode(to) )
                        return false ;
                    auto fNode = GetNode(from) ;
                    for( const auto  & eId : fNode.edge_ids)
                    {
                        if(  GetEdge(eId) == tmp )
                            return true ;
                    }
                    return false ;
                }

                size_t EdgesSize() const 
                {
                    return edges.size();
                }

                size_t NodesSize() const 
                {
                    return nodes.size();
                }

                void PrintAsDOT(std::ostream & out)
                {
                    out<<Edge::DOTHead()<<std::endl;
                    for( const auto & e : edges )
                    {
                        out<<"\t"<<e.ToString()<<std::endl;
                    }
                    out<<"}"<<std::endl;
                }
            };

        //
        //
        // Undirected Graph 
        //
        //
        template<class TNode 
            , class TEdge 
            , class TNodes = std::map<typename TNode::NodeNodeId , TNode>
            , class TEdges = std::vector<TEdge> 
            >
            struct ListGraph : public ListGraphBasic<TNode , TEdge, TNodes , TEdges>
            {
                typedef TNode Node ;
                typedef TEdge Edge ;
                typedef TNodes Nodes ;
                typedef TEdges Edges ;

                typedef typename Node::NodeNodeId NodeId ;
                typedef typename Edge::EdgeEdgeId EdgeId ;

                typedef ListGraphBasic<TNode , TEdge, TNodes , TEdges> Basic;

                // DO NOT use this if your TEdge has more member variable.
                void AddEdge( const typename Basic::NodeId & from ,const typename Basic::NodeId & to )
                {
                    // Make a new edge
                    TEdge tmp ;
                    tmp.from = from ;
                    tmp.to = to ;
                    AddEdge(tmp);
                }

                void AddEdge(const TEdge &tmp)
                {
                    auto from = tmp.from ;
                    auto to = tmp.to ;

                    if( Basic::CheckEdge(from , to ) )
                        return ;
                    if( ! Basic::HasNode( from ) )
                        Basic::AddNode(from) ;
                    if ( ! Basic::HasNode( to ) )
                        Basic::AddNode(to );

                    size_t nId = Basic::edges.size() ;
                    Basic::edges.push_back(tmp);
                    Basic::edges.rbegin()->id = nId;
                    Basic::GetNode(from).AddEdge(nId);
                    Basic::GetNode(to).AddEdge(nId);
                }

                template<class Me>
                Me SubGraph(const std::set<NodeId>& subs) const 
                {
                    Me ret ;

                    for(const auto & id : subs )
                    {
                        if( !Basic::HasNode( id ) )
                            continue;
                        const auto & node = Basic::GetNode(id);
                        ret.AddNode(node);
                        ret.GetNode(id).edge_ids.clear();
                    }

                    for(const auto & id : subs )
                    {
                        if( !Basic::HasNode( id ) )
                            continue;
                        const auto & node = Basic::GetNode(id);
                        for( const auto & edge_id : node.edge_ids)
                        {
                            const auto & edge = GetEdge(edge_id);
                            if( ret.HasNode( edge.from) && ret.HasNode(edge.to) )
                                ret.AddEdge(edge);
                        }
                    }
                    return ret;
                };
            };

        //
        //
        // Digraph
        //
        //
        template<class TNode 
            , class TEdge 
            , class TNodes = std::map<typename TNode::NodeNodeId , TNode>
            , class TEdges = std::vector<TEdge> 
            >
            struct ListDigraph  : public ListGraphBasic<TNode , TEdge, TNodes , TEdges>
            {
                typedef TNode Node ;
                typedef TEdge Edge ;
                typedef TNodes Nodes ;
                typedef TEdges Edges ;

                typedef typename Node::NodeNodeId NodeId ;
                typedef typename Edge::EdgeEdgeId EdgeId ;

                typedef ListGraphBasic<TNode , TEdge, TNodes , TEdges> Basic;

                // DO NOT use this if your TEdge has more member variable.
                void AddEdge( const typename Basic::NodeId & from , const typename Basic::NodeId & to )
                {
                    // Make a new edge
                    TEdge tmp ;
                    tmp.from = from ;
                    tmp.to = to ;
                    AddEdge(tmp);
                }

                void AddEdge(const TEdge &tmp)
                {
                    auto from = tmp.from ;
                    auto to = tmp.to ;

                    if( Basic::CheckEdge(from , to ) )
                        return ;
                    if( ! Basic::HasNode( from ) )
                        Basic::AddNode(from) ;
                    if ( ! Basic::HasNode( to ) )
                        Basic::AddNode(to );


                    size_t nId = Basic::edges.size() ;
                    Basic::edges.push_back(tmp);
                    Basic::edges.rbegin()->id = nId;
                    Basic::GetNode(from).AddEdge(nId);
                }

                template<class Me>
                Me SubGraph(const std::set<NodeId>& subs) const 
                {
                    Me ret ;

                    for(const auto & id : subs )
                    {
                        if( !Basic::HasNode( id ) )
                            continue;
                        const auto & node = Basic::GetNode(id);
                        ret.AddNode(node);
                        ret.GetNode(id).edge_ids.clear();
                    }

                    for(const auto & id : subs )
                    {
                        if( !Basic::HasNode( id ) )
                            continue;
                        const auto & node = Basic::GetNode(id);
                        for( const auto & edge_id : node.edge_ids)
                        {
                            const auto & edge = Basic::GetEdge(edge_id);
                            if( ret.HasNode( edge.from) && ret.HasNode(edge.to) )
                                ret.AddEdge(edge);
                        }
                    }
                    return ret;
                }

            };

        // Must be implement before use
        template<class BaseGraph 
            , class NodeId
            , class EdgeId
            , class Node1
            , class Edge1
            >
            struct GraphAccessBase
            {
                typedef NodeId                       GraphNodeId ;
                typedef EdgeId                       GraphEdgeId ;

                typedef Node1                        Node;
                typedef Edge1                        Edge;

                Node & AccessNode(const NodeId &) { assert(0) ; }
                Edge & AccessEdge(const EdgeId & , const GraphNodeId &){ assert(0); }

                BaseGraph * base ;

                std::map<GraphNodeId , Node > nodes ;
                std::map<GraphEdgeId , Edge > edges ;
                std::map<GraphNodeId , Edge > edges_ends;

                template<class Value>
                    void GetAttr(const Node & , const Edge & ,
                            const std::string & , Value & ) {assert(0);} 

            };

        template<class GraphAccess>
            struct EdgeIterator : public std::iterator<std::forward_iterator_tag,int>
        {
            public:

                typedef typename GraphAccess::GraphEdgeId Id;
                typedef typename GraphAccess::Edge        Edge;

                EdgeIterator() : curr(NULL) ,accessor( NULL ) { }

                EdgeIterator(const Edge & e , GraphAccess & acc )
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
                            curr = &(accessor->AccessEdge(next, curr->from));
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
                GraphAccess * accessor ;
        };


        template<class GraphAccess , class traits , class SearchNode  >
            struct PathEndHelperBase
            {
                typedef typename GraphAccess::GraphNodeId NodeId;
                typedef typename GraphAccess::GraphEdgeId EdgeId;
                typedef typename GraphAccess::Node        Node;
                typedef typename GraphAccess::Edge        Edge;
                typedef SearchNode                        SNode;
                typedef GraphAccess                       Access;
                typedef traits                            traisId;

                void Start() {assert(0);}
                void AddNode(const Node & , const SearchNode &) {assert(0);}
                void AddEdge(const Edge & ) {assert(0);}
                void PopEdge() {assert(0);}
                void PopNode() {assert(0);}
                bool IsEnd() const { assert(0) ;} ;

            };
    } // GRAPH
} // BGIQD 

#endif //__ALGORITHM_GRAPH_H__
