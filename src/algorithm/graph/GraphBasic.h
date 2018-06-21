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

                std::string ToString() const
                {
                    std::ostringstream ost;
                    ost<<from<<"\t--\t"<<to<<" [ id = "<<id << "]";
                    return ost.str();
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
                std::string ToString() const
                {
                    std::ostringstream ost;
                    ost<<BaseType::from<<"\t->\t"<<BaseType::to<<" [ id = "<<BaseType::id << "]";
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

                void AddEdge(const NodeEdgeId id )
                {
                    edge_ids.insert(id);
                }

                bool HasEdge(const NodeEdgeId id )
                {
                    edge_ids.find(id) != edge_ids.end() ; 
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

                void AddNode(NodeId id )
                {
                    nodes[id].id = id ;
                    //nodes[id].edge_id = Edge::invalid ;
                }

                void AddNode(const Node & n )
                {
                    nodes[n.id] = n ;
                }

                Node & GetNode( NodeId id )
                {
                    return nodes[id] ;
                }

                bool HasNode( NodeId id )
                {
                    return nodes.find(id) != nodes.end() ;
                }

                Edge & GetEdge( EdgeId id )
                {
                    return edges[id] ;
                }

                bool CheckEdge(NodeId from , NodeId to )
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

                void PrintAsDOT()
                {
                    std::cout<<Edge::DOTHead()<<std::endl;
                    for( const auto & e : edges )
                    {
                        std::cout<<"\t"<<e.ToString()<<std::endl;
                    }
                    std::cout<<"}"<<std::endl;
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
                void AddEdge( typename Basic::NodeId from , typename Basic::NodeId to )
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

                    Basic::edges.push_back(tmp);

                    size_t nId = Basic::edges.size() ;
                    Basic::edges.rbegin()->id = nId;
                    Basic::GetNode(from).AddEdge(nId);
                    Basic::GetNode(to).AddEdge(nId);
                }
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
                void AddEdge( typename Basic::NodeId from , typename Basic::NodeId to )
                {
                    if( Basic::CheckEdge(from , to ) )
                        return ;
                    if( ! Basic::HasNode( from ) )
                        Basic::AddNode(from) ;
                    if ( ! Basic::HasNode( to ) )
                        Basic::AddNode(to );
                    size_t nId = Basic::edges.size() ;
                    TEdge tmp ;
                    tmp.id = nId ;
                    tmp.from = from ;
                    tmp.to = to ;
                    // Make a new edge
                    Basic::edges.push_back(tmp);
                    Basic::GetNode(from).AddEdge(nId);
                }
            };

    } // GRAPH
} // BGIQD 

#endif //__ALGORITHM_GRAPH_H__
