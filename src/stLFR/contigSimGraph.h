#ifndef __STLFR_CONTIGSIMGRAPH_H__
#define __STLFR_CONTIGSIMGRAPH_H__

#include <vector>
#include <set>
#include <map>

#include "algorithm/graph/GraphBasic.h"
#include "algorithm/graph/MinTree.h"
#include "algorithm/graph/GraphTrunk.h"
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
            typedef BGIQD::GRAPH::ListGraph<Node,Edge> Basic ;

            void AddEdgeSim( unsigned int from , unsigned int to , float sim)
            {
                Edge tmp ;
                tmp.from = from ;
                tmp.to = to ;
                tmp.sim = sim ;
                Basic::AddEdge(tmp);
            }

            typedef BGIQD::GRAPH::MinTreeHelper<ContigSimGraph, float , EdgeAttr> MTHelper;
            typedef BGIQD::GRAPH::TrunkHelper< ContigSimGraph> TKHelper;

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

            std::vector<Basic::NodeId> TrunkLinear(const ContigSimGraph & mintree)
            {
                return TKHelper::LinearTrunk( mintree );
            }
        };
    }
}

#endif //__STLFR_CONTIGSIMGRAPH_H__
