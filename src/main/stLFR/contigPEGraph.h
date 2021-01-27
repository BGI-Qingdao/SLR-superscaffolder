#ifndef __STLFR_CONTIGPEGRAPH_H__
#define __STLFR_CONTIGPEGRAPH_H__

#include "utils/graph/GraphBasic.h"
#include "utils/flags/flags.h"

#include <sstream>

namespace BGIQD {
    namespace stLFR {

        struct ContigNode : public BGIQD::GRAPH::IGraphNodeBasic<unsigned int , long >
        {
            typedef BGIQD::GRAPH::IGraphNodeBasic<unsigned int , long >  Basic;
            int contigLen ;
            long edge_id;
        };

        struct PEEdge : public BGIQD::GRAPH::IDigraphEdgeBase<unsigned int ,long>
        {
            int count ;
            int len;
            long next ; // to use the SPFsearch code .
            //FLAGS_INT 
            //ADD_A_FLAG( 0 , R2R );

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

        struct ContigPEGraph : public BGIQD::GRAPH::ListDigraph<ContigNode , PEEdge >
        {
            typedef BGIQD::GRAPH::ListDigraph<ContigNode , PEEdge > Basic;
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

            void MakeEdgeNext(const NodeId & id)
            {
                auto & node = GetNode(id);
                node.edge_id = Edge::invalid ;
                if( node.edge_ids.empty() )
                {
                    return ;
                }

                std::vector<std::tuple<int,long> >  edges;

                for( auto i : node.edge_ids )
                {
                    auto & edge = GetEdge(i);
                    edges.push_back(std::make_tuple(edge.count, i) );
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
