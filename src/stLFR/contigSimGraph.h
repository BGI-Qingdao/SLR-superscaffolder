#ifndef __STLFR_CONTIGSIMGRAPH_H__
#define __STLFR_CONTIGSIMGRAPH_H__

#include <vector>
#include <set>
#include <map>

#include "algorithm/graph/GraphBasic.h"
#include "algorithm/graph/MinTree.h"
#include "algorithm/graph/GraphTrunk.h"
#include "algorithm/disjoin_set/disjoin_set.h"

namespace BGIQD {
    namespace stLFR {

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
                ost<<" id= "<<id<<" sim= "<<sim;
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
            typedef BGIQD::GRAPH::TrunkHelper< ContigSimGraph> TKHelper;

            typedef BGIQD::Algorithm::DisJoin_Set<NodeId> DJ_Sets;

            ContigSimGraph  MinTree() const 
            {
                EdgeAttr attr ;
                MTHelper helper;
                return helper.MinTree(*this , attr);
            };

            ContigSimGraph TrunkFromMinTree(const ContigSimGraph & mintree)
            {
                return TKHelper::Trunk(mintree);
            }

            static std::map<NodeId , ContigSimGraph>  UnicomGraph(const ContigSimGraph & mintree)
            {
                std::map<NodeId , ContigSimGraph> ret;
                DJ_Sets dj_sets;
                for( const auto & edge : mintree.edges )
                {
                    dj_sets.AddConnect(edge.from , edge.to);
                }
                for( const auto & edge : mintree.edges )
                {
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
            std::vector<Basic::NodeId> TrunkLinear(const ContigSimGraph & mintree)
            {
                return TKHelper::LinearTrunk( mintree );
            }
        };
    }
}

#endif //__STLFR_CONTIGSIMGRAPH_H__
