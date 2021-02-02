#ifndef __STLFR_CONTIGSIMGRAPH_H__
#define __STLFR_CONTIGSIMGRAPH_H__

#include <vector>
#include <set>
#include <map>

#include "utils/graph/GraphBasic.h"
#include "utils/graph/mst/MinTree.h"
#include "utils/graph/DisJointSet.h"
#include "utils/graph/GraphTipRemove.h"

namespace BGIQD {
    namespace stLFR {

        template<class TEdge>
            struct TrunkNode
            {
                typedef typename TEdge::EdgeNodeId NodeId;
                typedef typename TEdge::EdgeEdgeId EdgeId;
                std::set<EdgeId> edges;

                int level;
                NodeId prev;
                NodeId prev_2;
                NodeId id;
                bool marked;
            };

        template<class TListGraph ,
            class TNode = TrunkNode<typename TListGraph::Edge> 
                >
                struct TrunkHelper
                {
                    typedef TListGraph ListGraph ;
                    typedef TNode TrunkNode ;

                    typedef typename TListGraph::Node Node;
                    typedef typename TListGraph::Edge Edge;
                    typedef typename Edge::EdgeEdgeId EdgeId;
                    typedef typename Edge::EdgeNodeId NodeId;


                    // Make sure input is a trunk
                    static std::vector<NodeId> LinearTrunk(const ListGraph &  base )
                    {
                        std::vector<NodeId> ret ;
                        NodeId starter ;
                        for( const auto & pair : base.nodes )
                        {
                            auto & node = pair.second ;
                            if( node.edge_ids.size() == 1 )
                            {
                                starter = node.id ;
                                break;
                            }
                            else
                            {
                                assert( node.edge_ids.size() == 2 );
                            }
                        }
                        ret.push_back( starter ) ;

                        auto & node = base.GetNode(starter);;
                        auto edge_id = *(node.edge_ids.begin());
                        auto & edge = base.GetEdge(edge_id);
                        NodeId next = edge.OppoNode(starter);
                        NodeId curr = starter;
                        while(1)
                        {
                            auto & next_node = base.GetNode(next) ;
                            ret.push_back(next);
                            if( next_node.edge_ids.size() == 1 )
                            {
                                break;
                            }
                            else
                            {
                                assert( next_node.edge_ids.size() == 2 );
                                auto edge_id1 = *(next_node.edge_ids.begin());
                                auto edge_id2 = *(std::next(next_node.edge_ids.begin()));
                                auto &edge1 = base.GetEdge(edge_id1);
                                auto &edge2 = base.GetEdge(edge_id2);
                                if( edge1.OppoNode(next) != curr)
                                {
                                    curr = next ;
                                    next = edge1.OppoNode(next) ;
                                }
                                else
                                {
                                    assert(edge2.OppoNode(next) != curr);
                                    curr = next ;
                                    next = edge2.OppoNode(next);
                                }
                            }
                        }
                        return ret ;
                    }
                };

        struct Node : public BGIQD::GRAPH::IGraphNodeBasic<unsigned int , long> 
        {
            enum Type 
            {
                Unknow = 0 ,
                Tip = 1 ,
                Single = 2 ,
                Linear = 3 ,
                Junction = 4  ,
            };

            Type type ;
        };

        struct Edge : public BGIQD::GRAPH::IGraphEdgeBasic<unsigned int , long >
        {
            float sim ;

            std::string AttrString() const
            {
                std::ostringstream ost;
                //ost<<" id= "<<id<<" sim= "<<sim;
                ost<<"label=\""<<sim<<"\"";
                return ost.str();
            }
            std::string ToString() const
            {
                std::ostringstream ost;
                ost<<from<<"\t--\t"<<to<<" [ "<<AttrString()<<" ]";
                return ost.str();
            }
        };

        struct EdgeAttr
        {
            float GetValue(const Edge & e ) const 
            {
                return 1.0f - e.sim ;
            }
        };

        struct ContigSimGraph : public BGIQD::GRAPH::ListGraph<Node,Edge>
        {

            bool use_salas ;

            typedef BGIQD::GRAPH::ListGraph<Node,Edge> Basic ;
            typedef Basic::NodeId NodeId;
            typedef Basic::EdgeId EdgeId;

            void AddEdgeSim( unsigned int from , unsigned int to , float sim)
            {
                // Force make linear in salas strategy.
                if( use_salas )
                    if( ( Basic::HasNode(from) && Basic::GetNode(from).edge_ids.size() >= 2) 
                            || ( Basic::HasNode(to) &&Basic::GetNode(to).edge_ids.size()>=2) )
                        return ;
                Edge tmp ;
                tmp.from = from ;
                tmp.to = to ;
                tmp.sim = sim ;
                Basic::AddEdge(tmp);
            }

            typedef BGIQD::GRAPH::MinTreeHelper<ContigSimGraph, float , EdgeAttr> MTHelper;
            typedef BGIQD::GRAPH::TipRemoveHelper<ContigSimGraph> TipHelper ;
            typedef TipHelper::TipRemoveResult TipRemoveResult;

            typedef BGIQD::Algorithm::DisJoin_Set<NodeId> DJ_Sets;
            typedef TrunkHelper< ContigSimGraph> TKHelper;


            static TipRemoveResult RemoveTip_n2( ContigSimGraph & mintree )
            {
                TipHelper tip_helper ;

                auto tip_checker = [](const TipHelper::tip & t ) -> bool
                {
                    return t.size() < 2 ;
                };
                tip_helper.Init(tip_checker);
                return tip_helper.DeepTipRemove(mintree);
            }


            static std::vector<NodeId> DetectJunctions( const ContigSimGraph & mintree)
            {
                std::vector<NodeId> ret ;
                for( const auto & pair : mintree.nodes )
                {
                    const auto & node = pair.second ;
                    if( node.EdgeNum() > 2 )
                    {
                        ret.push_back( node.id ) ;
                    }
                }
                return ret;
            }

            ContigSimGraph  MinTree() const 
            {
                EdgeAttr attr ;
                MTHelper helper;
                return helper.MinTree(*this , attr);
            };

            static std::map<NodeId , ContigSimGraph>  UnicomGraph(const ContigSimGraph & mintree)
            {
                std::map<NodeId , ContigSimGraph> ret;
                DJ_Sets dj_sets;
                for( const auto & edge : mintree.edges )
                {
                    if( edge.IsValid())
                        dj_sets.AddConnect(edge.from , edge.to);
                }
                for( const auto & edge : mintree.edges )
                {
                    if( ! edge.IsValid() )
                        continue;
                    auto rep = dj_sets.GetGroup(edge.from);
                    if( ! ret[rep].HasNode(edge.from) )
                    {
                        ret[rep].AddNode(mintree.GetNode(edge.from));
                        ret[rep].GetNode(edge.from).edge_ids.clear();
                    }
                    if( ! ret[rep].HasNode(edge.to) )
                    {
                        ret[rep].AddNode(mintree.GetNode(edge.to));
                        ret[rep].GetNode(edge.to).edge_ids.clear();
                    }
                    ret[rep].AddEdge(edge);
                }
                return ret ;
            }
            static std::vector<Basic::NodeId> TrunkLinear(const ContigSimGraph & mintree)
            {
                return TKHelper::LinearTrunk( mintree );
            }
            void PrintAsDOT(std::ostream & out) const
            {
                out<<Edge::DOTHead()<<std::endl;
                std::set<int> multis;
                for( const auto & e : edges )
                {
                    out<<"\t"<<e.ToString()<<std::endl;
                    if( GetNode(e.from).EdgeNum() > 2 )
                        multis.insert(e.from);
                    if(  GetNode(e.to).EdgeNum()>2 )
                        multis.insert(e.to);
                }
                for( int i : multis )
                {
                    out<<i<< " [ shape=box]  "<<std::endl;
                }
                out<<"}"<<std::endl;
            }
        };
    }
}

#endif //__STLFR_CONTIGSIMGRAPH_H__
