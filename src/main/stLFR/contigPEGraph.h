#ifndef __STLFR_CONTIGPEGRAPH_H__
#define __STLFR_CONTIGPEGRAPH_H__

#include "utils/graph/GraphBasic.h"
#include "utils/misc/flags.h"

#include <sstream>

/**********************************************************
 *
 * @Brief  :
 *
 *     The contig-PE graph used by
 *          + PEGraph
 *          + FillTrunkByPE
 *
 * *******************************************************/

namespace BGIQD {
    namespace stLFR {

        // Node of PE graph.
        // To achive a new edge iterator , we add edge_id to store the strongest PE link edge
        struct ContigNode : public BGIQD::GRAPH::IGraphNodeBasic<unsigned int , long >
        {
            typedef BGIQD::GRAPH::IGraphNodeBasic<unsigned int , long >  Basic;
            int contigLen ;
            long edge_id;
        };

        // Edge of PE graph.
        struct PEEdge : public BGIQD::GRAPH::IDigraphEdgeBase<unsigned int ,long>
        {
            int count ; // PE linkage count.
            int len   ; // gap size estimated by PE.
            long next ; // to use the Search code . This search code will be replaced soon.

            std::string AttrString() const
            {
                std::ostringstream ost;
                ost<<" count= "<<count<<" len= "<<len;//<<" flags= "<<flags;
                return ost.str();
            }
            std::string ToString() const
            {
                std::ostringstream ost;
                ost<<from<<"\t->\t"<<to<<" [ "<<AttrString()<<" ]";
                return ost.str();
            }

            void InitFromString(const std::string & line )
            {
                //\t317165\t->\t346282 [  count= 11 len= -153 ]
                char tmp_char ;
                std::string tmp_str;
                std::istringstream ist(line);
                ist>>from>>tmp_str>>to>>tmp_char>>tmp_str>>count>>tmp_str>>len;
            }
        };

        // contig-PE graph for local scaffolding.
        struct ContigPEGraph : public BGIQD::GRAPH::Digraph<ContigNode , PEEdge >
        {
            typedef BGIQD::GRAPH::Digraph<ContigNode , PEEdge > Basic;
            // Must guarantee contigId +1 is the compelete reverse contig if this one.
            void AddNode( unsigned int contigId , int len)
            {
                Node tmp ;
                tmp.contigLen = len ;
                tmp.id = contigId ;
                Basic::AddNode(tmp);
                tmp.id = contigId +1 ;
                Basic::AddNode(tmp);
            }

            // Sort edges of one node by the PE linkage count and chained in list.
            void MakeEdgeNext(const NodeId & id)
            {
                auto & node = GetNode(id);
                node.edge_id = Edge::invalid ;
                if( node.EdgeNum() == 0)
                {
                    return ;
                }

                std::vector<std::tuple<int,long> >  edges;
                ContigPEGraph::Node::NodeEdgeIdIterator begin, end;
                std::tie(begin,end) = node.GetEdges();
                for(auto x = begin ; x!= end ; x++ )
                {
                    auto & edge = GetEdge(*x);
                    edges.push_back(std::make_tuple(edge.count, *x) );
                }

                std::sort(edges.rbegin() ,edges.rend());
                node.edge_id = std::get<1>(*edges.begin());

                EdgeId next = Edge::invalid ;
                for( auto i = edges.rbegin() ; i != edges.rend() ; i =  std::next(i) )
                {
                    EdgeId eid = std::get<1>(*i);
                    auto & edge = GetEdge(eid);
                    if( next != Edge::invalid )
                    {
                        edge.next = next ;
                    }
                    else
                    {
                        edge.next = Edge::invalid ;
                    }
                    next = eid ;
                }
            }

            void AddEdge( unsigned int from , unsigned int to ,int len, int count ) 
            {
                if( ! ( Basic::HasNode( from ) && Basic::HasNode(to) ) )
                    return ;
                Edge edge;
                edge.from = from ;
                edge.to = to ; 
                edge.count = count ;
                edge.len = len ;
                edge.next = Edge::invalid ;
                Basic::AddEdge(edge);
            }

            ContigPEGraph SubGraph(const std::set<NodeId>& subs) const
            {
                return Basic::SubGraph<ContigPEGraph>(subs);
            }
            void AddEdge(const Edge & edge)
            {
                Basic::AddEdge(edge);
            }

            void AddNode(const Node & tmp)
            {
                Basic::AddNode(tmp);
            }
        };
    }
}
#endif
