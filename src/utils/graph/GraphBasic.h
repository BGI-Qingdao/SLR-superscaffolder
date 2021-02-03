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
#include <tuple>
/**********************************************************
 *
 * @Brief
 *  A simple implement of Graph structure.
 *      + one un-directed graph
 *      + one directed graph
 *
 * *******************************************************/

namespace BGIQD {
    namespace GRAPH {

        // An un-directed edge
        template<class NodeId ,  class EdgeId >
            struct IGraphEdgeBasic
            {
                typedef NodeId EdgeNodeId;
                typedef EdgeId EdgeEdgeId;

                EdgeId      id ;

                NodeId      from ;
                NodeId      to ;

                bool operator == (const  IGraphEdgeBasic & i ) const
                {
                    return ( (from == i.from && to == i.to ) 
                            || ( from == i.to && to == i.from ) ) ;
                }

                bool operator != ( const  IGraphEdgeBasic & i ) const
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
                // This may be used for lazy-deletion
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

        // A directed edge
        template<class NodeId ,  class EdgeId >
            struct IDigraphEdgeBase : public IGraphEdgeBasic<NodeId, EdgeId>
            {
                typedef IGraphEdgeBasic<NodeId, EdgeId> BaseType ;

                typedef NodeId EdgeNodeId;
                typedef EdgeId EdgeEdgeId;

                bool operator == (const  IDigraphEdgeBase & i ) const
                {
                    return ( (BaseType::from == i.from && BaseType::to == i.to ));
                }

                bool operator != ( const IDigraphEdgeBase & i ) const
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

        // A node struct.
        //
        // This node contain all ids of its out-edges but struct edge is not stored by node.
        //
        // Notice : for undirected graph, all edges are out-edge.
        template<class NodeId ,  class EdgeId >
            struct IGraphNodeBasic
            {
                typedef NodeId NodeNodeId;
                typedef EdgeId NodeEdgeId;

                NodeId                         id ;

                typedef typename std::set<EdgeId>::iterator NodeEdgeIdIterator;

                bool HasEdge(const NodeEdgeId &id )
                {
                    return edge_ids.find(id) != edge_ids.end() ; 
                }

                size_t EdgeNum() const 
                {
                    return edge_ids.size() ;
                }

                std::pair<NodeEdgeIdIterator,NodeEdgeIdIterator> GetEdges() const
                {
                    return std::make_pair(edge_ids.begin() , edge_ids.end());
                }

                void AddEdge(const NodeEdgeId &id )
                {
                    edge_ids.insert(id);
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

                protected:
                std::set<EdgeId>               edge_ids;
            };

        //
        //
        // GraphBaisc.
        //
        // This graph can be iteratered by both nodes and edges.
        //
        template<class TNode 
            , class TEdge
            , class TNodes = std::map<typename TNode::NodeNodeId , TNode>
            , class TEdges = std::vector<TEdge> 
            >
            struct IGraphBasic
            {
                typedef TNode Node ;
                typedef TEdge Edge ;
                typedef TNodes Nodes ;
                typedef TEdges Edges ;

                typedef typename Node::NodeNodeId NodeId ;
                typedef typename Edge::EdgeEdgeId EdgeId ;

                Nodes nodes ;
                Edges edges ;

                // Access node :
                Node & GetNode(const NodeId & id )
                {
                    return nodes[id] ;
                }

                const Node & GetNode(const NodeId & id ) const 
                {
                    return nodes.at(id) ;
                }

                bool HasNode( const NodeId & id ) const
                {
                    return nodes.find(id) != nodes.end() ;
                }

                // Access Edge:
                //
                size_t NodesSize() const 
                {
                    return nodes.size();
                }

                size_t EdgesSize() const 
                {
                    return edges.size();
                }

                Edge & GetEdge( const EdgeId & id )
                {
                    return edges[id] ;
                }

                const Edge & GetEdge( const EdgeId &id ) const
                {
                    return edges[id] ;
                }

                bool CheckEdge(const NodeId & from ,const  NodeId & to ) const
                {
                    Edge tmp ;
                    tmp.from = from ;
                    tmp.to = to ;
                    if( ! HasNode( from ) || !HasNode(to) )
                        return false ;
                    auto fNode = GetNode(from) ;
                    typename Node::NodeEdgeIdIterator begin, end ;
                    std::tie(begin,end) = fNode.GetEdges();
                    for( auto x = begin ; x != end ; x++)
                    {
                        if(  GetEdge(*x) == tmp )
                            return true ;
                    }
                    return false ;
                }
                // Modify node
                void AddNode(const Node & n )
                {
                    nodes[n.id] = n ;
                }

                bool RemoveNode( const NodeId & id )
                {
                    if ( nodes.find(id) != nodes.end() )
                    {
                        auto & n1 = GetNode( id );
                        typename Node::NodeEdgeIdIterator begin, end ;
                        std::tie(begin,end) = n1.GetEdges();
                        for( auto x = begin ; x != end ; x++)
                        {
                            RemoveEdge(*x);
                        }
                        nodes.erase(id);
                        return true ;
                    }
                    else
                        return false ;
                }

                // Modify edge
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

                // Print graph in DOT format
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
            struct Graph : public IGraphBasic<TNode , TEdge, TNodes , TEdges>
            {
                typedef TNode Node ;
                typedef TEdge Edge ;
                typedef TNodes Nodes ;
                typedef TEdges Edges ;

                typedef typename Node::NodeNodeId NodeId ;
                typedef typename Edge::EdgeEdgeId EdgeId ;

                typedef IGraphBasic<TNode , TEdge, TNodes , TEdges> Basic;

                // Add edge and manager correspondence nodes
                void AddEdge(const TEdge &tmp)
                {
                    auto from = tmp.from ;
                    auto to = tmp.to ;

                    if( Basic::CheckEdge(from , to ) )
                        return ;
                    assert( Basic::HasNode( from ) == true);
                    assert( Basic::HasNode( to ) == true ) ;

                    size_t nId = Basic::edges.size() ;
                    Basic::edges.push_back(tmp);
                    Basic::edges.rbegin()->id = nId;
                    Basic::GetNode(from).AddEdge(nId);
                    Basic::GetNode(to).AddEdge(nId);
                }
                // get a subgraph based on given node set
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
                        ret.GetNode(id).CleanEdges();
                    }

                    for(const auto & id : subs )
                    {
                        if( !Basic::HasNode( id ) )
                            continue;
                        const auto & node = Basic::GetNode(id);
                        typename Node::NodeEdgeIdIterator begin, end ;
                        std::tie(begin,end) = node.GetEdges();
                        for(auto i = begin; i!=end; i++)
                        {
                            const auto & edge = GetEdge(*i);
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
            struct Digraph  : public IGraphBasic<TNode , TEdge, TNodes , TEdges>
            {
                typedef TNode Node ;
                typedef TEdge Edge ;
                typedef TNodes Nodes ;
                typedef TEdges Edges ;

                typedef typename Node::NodeNodeId NodeId ;
                typedef typename Edge::EdgeEdgeId EdgeId ;

                typedef IGraphBasic<TNode , TEdge, TNodes , TEdges> Basic;

                // Add edge and manager correspondence nodes
                void AddEdge(const TEdge &tmp)
                {
                    auto from = tmp.from ;
                    auto to = tmp.to ;

                    if( Basic::CheckEdge(from , to ) )
                        return ;
                    assert( Basic::HasNode( from ) == true);
                    assert( Basic::HasNode( to ) == true ) ;

                    size_t nId = Basic::edges.size() ;
                    Basic::edges.push_back(tmp);
                    Basic::edges.rbegin()->id = nId;
                    Basic::GetNode(from).AddEdge(nId);
                }

                // get a subgraph based on given node set
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
                        ret.GetNode(id).CleanEdges();
                    }

                    for(const auto & id : subs )
                    {
                        if( !Basic::HasNode( id ) )
                            continue;
                        const auto & node = Basic::GetNode(id);
                        typename Node::NodeEdgeIdIterator begin, end ;
                        std::tie(begin,end) = node.GetEdges();
                        for(auto i = begin; i!=end; i++)
                        {
                            const auto & edge = Basic::GetEdge(*i);
                            if( ret.HasNode( edge.from) && ret.HasNode(edge.to) )
                                ret.AddEdge(edge);
                        }
                    }
                    return ret;
                }

            };

    } // GRAPH
} // BGIQD 

#endif //__ALGORITHM_GRAPH_H__
